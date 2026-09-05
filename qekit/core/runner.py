# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Ejecución por lotes de cálculos de pw.x.

Los módulos de barrido (convergencia, ecuación de estado, constantes
elásticas) generan decenas de cálculos pequeños. Este módulo los corre uno
tras otro, muestra el avance y recolecta los resultados.

Dos decisiones de diseño que importan en la práctica:

- **Se puede reanudar.** Antes de lanzar un cálculo se comprueba si ya
  terminó de verdad: no basta con que ponga `JOB DONE.`, tiene que tener un
  XML legible y haber convergido. Un cálculo que terminó sin converger se
  vuelve a lanzar en vez de darlo por bueno, que es como se cuela un punto
  malo en una curva sin que nadie lo note.
- **Un fallo no aborta el barrido.** Si un punto no converge, se anota y se
  sigue con los demás: es preferible obtener nueve puntos de diez y saber
  cuál falló, que perder toda la serie.
- **Se pueden correr varios a la vez.** pw.x escala mal: en una máquina de
  ocho hilos, cuatro cálculos de dos procesos terminan un barrido bastante
  antes que uno de ocho procesos repetido cuatro veces. Los puntos de un
  barrido son independientes y viven en carpetas distintas, así que la
  única precaución real es no pedir más procesos de los que hay.
"""

import os
import shutil
import subprocess
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path

from qekit import config as qcfg
from qekit.core import plataforma
from qekit.core import qeout

DONE_MARK = "JOB DONE."


def nucleos() -> int:
    """Hilos disponibles de verdad, respetando el cgroup si lo hay."""
    try:
        return len(os.sched_getaffinity(0))
    except AttributeError:                                  # pragma: no cover
        return os.cpu_count() or 1


def memoria_libre_gb():
    """GB de memoria disponible, en los tres sistemas.

    Antes solo leía /proc/meminfo y devolvía None en macOS y en Windows, así
    que el aviso de «vas a sobresuscribir la memoria» —justo el que evita
    que un barrido en paralelo tumbe el portátil— no aparecía nunca fuera de
    Linux. Cada sistema tiene su forma de contarlo y ninguna necesita
    instalar nada.
    """
    from qekit.core import plataforma
    if plataforma.LINUX:
        try:
            for linea in Path("/proc/meminfo").read_text().splitlines():
                if linea.startswith("MemAvailable:"):
                    return float(linea.split()[1]) / 1048576.0
        except OSError:                                     # pragma: no cover
            pass
        return None
    if plataforma.MACOS:                                    # pragma: no cover
        try:
            import subprocess as _sp
            salida = _sp.run(["vm_stat"], capture_output=True, text=True,
                             timeout=5).stdout
            tam = 4096
            libres = 0
            for linea in salida.splitlines():
                if "page size of" in linea:
                    tam = int(linea.split("page size of")[1].split()[0])
                for etiqueta in ("Pages free:", "Pages inactive:",
                                 "Pages speculative:"):
                    if linea.startswith(etiqueta):
                        libres += int(linea.split(":")[1].strip().rstrip("."))
            return libres * tam / 1073741824.0 if libres else None
        except Exception:                                   # noqa: BLE001
            return None
    if plataforma.WINDOWS:                                  # pragma: no cover
        try:
            import ctypes

            class _Mem(ctypes.Structure):
                _fields_ = [("dwLength", ctypes.c_ulong),
                            ("dwMemoryLoad", ctypes.c_ulong),
                            ("ullTotalPhys", ctypes.c_ulonglong),
                            ("ullAvailPhys", ctypes.c_ulonglong),
                            ("ullTotalPageFile", ctypes.c_ulonglong),
                            ("ullAvailPageFile", ctypes.c_ulonglong),
                            ("ullTotalVirtual", ctypes.c_ulonglong),
                            ("ullAvailVirtual", ctypes.c_ulonglong),
                            ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]

            m = _Mem()
            m.dwLength = ctypes.sizeof(_Mem)
            ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m))
            return m.ullAvailPhys / 1073741824.0
        except Exception:                                   # noqa: BLE001
            return None
    return None


def reparto(paralelo: int, nproc: int = None) -> tuple:
    """Cuántos cálculos a la vez y cuántos procesos MPI a cada uno.

    Si el usuario no fija --nproc, se reparten los hilos de la máquina entre
    los cálculos simultáneos. Poner ocho procesos en cada uno de cuatro
    cálculos en una máquina de ocho hilos no va cuatro veces más rápido: va
    más lento que uno solo, porque los treinta y dos procesos se pelean.
    """
    cfg = qcfg.load()
    hilos = nucleos()
    par = max(1, int(paralelo or 1))
    if nproc is not None:
        n = max(1, int(nproc))
    elif par > 1:
        n = max(1, hilos // par)
    else:
        n = max(1, int(cfg.get("nproc", 1) or 1))
    aviso = ""
    if par * n > hilos:
        aviso = (f"you asked for {par} simultaneous calculations with {n} process"
                 f"{'es' if n > 1 else ''} each = {par * n} processes, and the "
                 f"machine has {hilos} thread{'s' if hilos > 1 else ''}. "
                 f"Oversubscribing makes ALL of them slower.")
    return par, n, aviso


@dataclass
class Job:
    """Un cálculo de pw.x dentro de su propia carpeta."""

    name: str                       # etiqueta legible ("ecutwfc = 60 Ry")
    directory: Path
    input_file: str = "pw.in"
    output_file: str = "pw.out"
    meta: dict = field(default_factory=dict)   # datos del barrido (volumen, etc.)

    @property
    def input_path(self) -> Path:
        return Path(self.directory) / self.input_file

    @property
    def output_path(self) -> Path:
        return Path(self.directory) / self.output_file

    def is_done(self, estricto: bool = True) -> bool:
        """¿Terminó bien este cálculo?

        Con `estricto` (lo normal) no basta con que pw.x haya escrito
        "JOB DONE.": además tiene que existir un XML legible y el cálculo
        tiene que haber CONVERGIDO. Un scf que agota electron_maxstep sin
        converger termina limpiamente y escribe JOB DONE; si se da por
        bueno, al reanudar un barrido ese punto malo se queda para siempre
        y contamina la curva sin que nada avise.
        """
        try:
            if DONE_MARK not in self.output_path.read_text(errors="ignore"):
                return False
        except OSError:
            return False
        if not estricto:
            return True
        try:
            res = qeout.read_xml(str(self.directory))
        except Exception:                                   # noqa: BLE001
            return False
        # None significa que el XML no permitió demostrar convergencia. En
        # modo estricto solo una confirmación explícita cuenta como hecho.
        return res.converged is True


@dataclass
class JobResult:
    job: Job
    ok: bool = False
    seconds: float = 0.0
    skipped: bool = False
    error: str = ""
    result: "qeout.QEResult" = None      # None si falló

    @property
    def energy(self):
        return self.result.total_energy if self.result else None


def build_command(pw_cmd: str = None, nproc: int = None) -> list:
    """Comando para lanzar pw.x, en serie o con MPI.

    Se toma de la configuración (`pw_cmd`, `nproc`) y se puede sobrescribir.
    Si `pw_cmd` ya trae un lanzador (mpirun, srun...), se respeta tal cual.

    En Windows el binario se llama `pw.exe`, y el lanzador MPI que hay es
    `mpiexec` (MS-MPI), no `mpirun`. Las dos cosas se resuelven mirando qué
    existe de verdad en la máquina, no suponiendo Linux.
    """
    cfg = qcfg.load()
    cmd = (pw_cmd or cfg.get("pw_cmd") or "").strip()
    if not cmd:
        cmd = plataforma.nombres_ejecutable("pw")[0]
    n = int(nproc if nproc is not None else cfg.get("nproc", 1) or 1)

    if any(tok in cmd for tok in ("mpirun", "mpiexec", "srun")):
        return cmd.split()          # el usuario ya especificó cómo lanzarlo
    if n > 1:
        launcher = (cfg.get("mpi_cmd") or plataforma.lanzador_mpi())
        if launcher:
            return launcher.format(n=n).split() + cmd.split()
        # sin MPI en la máquina: correr en serie es mejor que fallar
    return cmd.split()


def resolver_ejecutable(exe: str) -> str:
    """La ruta real de un binario, probando .x y .exe y los sitios de siempre.

    `shutil.which` solo mira el PATH. En macOS con Homebrew y en muchas
    instalaciones a mano de Quantum ESPRESSO los binarios no están en el
    PATH, y decirle a alguien «no encuentro pw.x» en una máquina donde está
    instalado es la peor forma de recibirle.
    """
    directo = shutil.which(exe)
    if directo:
        return directo
    p = Path(exe)
    if p.is_absolute() or len(p.parts) > 1:
        for nombre in plataforma.nombres_ejecutable(p.name):
            cand = p.with_name(nombre)
            if cand.exists():
                return str(cand)
        return None
    return plataforma.buscar_ejecutable(p.name,
                                        plataforma.dirs_probables_qe())


def check_available(pw_cmd: str = None, nproc: int = None) -> str:
    """Comprueba que el ejecutable exista; devuelve la ruta encontrada."""
    cmd = build_command(pw_cmd, nproc)
    exe = cmd[_qe_executable_index(cmd)]
    found = resolver_ejecutable(exe)
    if not found:
        sug = ""
        if plataforma.WINDOWS:
            sug = ("\nOn Windows the Quantum ESPRESSO binaries are called "
                   "pw.exe, not pw.x.\nIf you have it in WSL, run Olla-DFT "
                   "inside WSL: there it will see the Linux pw.x.")
        raise FileNotFoundError(
            f"executable '{exe}' not found.\n"
            "Install Quantum ESPRESSO or set the path with:\n"
            f"  olla-dft config set pw_cmd /path/to/{exe}\n"
            "You can also generate the inputs without --run and run them yourself."
            + sug
        )
    return found


def _qe_executable_index(cmd: list) -> int:
    """Índice de pw.x/pw.exe en un comando, también detrás de MPI/Slurm."""
    if not cmd:
        raise ValueError("the Quantum ESPRESSO command is empty")
    if Path(cmd[0]).name.lower() not in ("mpirun", "mpiexec", "srun"):
        return 0
    for i in range(1, len(cmd)):
        nombre = Path(cmd[i]).name.lower()
        if nombre in ("pw", "pw.x", "pw.exe"):
            return i
    # build_command añade el comando de pw al final. Este respaldo conserva
    # ejecutables con un nombre personalizado, por ejemplo pw-custom.
    return len(cmd) - 1


def run_one(job: Job, cmd: list, timeout: float = None,
            rehacer: bool = False) -> JobResult:
    """Lanza un cálculo y devuelve su resultado."""
    res = JobResult(job=job)
    if not rehacer and job.is_done():
        res.ok, res.skipped = True, True
        try:
            res.result = qeout.read_xml(str(job.directory))
        except Exception:
            res.result = None
        return res

    start = time.time()
    try:
        with open(job.input_path) as fin, open(job.output_path, "w") as fout:
            proc = subprocess.run(cmd, stdin=fin, stdout=fout,
                                  stderr=subprocess.STDOUT,
                                  cwd=str(job.directory), timeout=timeout)
        res.seconds = time.time() - start
        if proc.returncode != 0:
            res.error = f"pw.x exited with code {proc.returncode}"
            # Un código de salida a secas no dice nada. Si la salida trae
            # una causa reconocible, se pega aquí para no tener que abrir
            # el archivo — sobre todo en un barrido de 25 cálculos.
            try:
                pista = failure_hint(job.output_path.read_text(errors="ignore"))
            except OSError:
                pista = ""
            if pista:
                res.error += f" — {pista}"
            return res
    except subprocess.TimeoutExpired:
        res.seconds = time.time() - start
        res.error = "the time limit was exceeded"
        return res
    except OSError as exc:
        res.error = str(exc)
        return res

    if not job.is_done(estricto=False):
        res.error = "the calculation did not finish (check the output)"
        try:
            pista = failure_hint(job.output_path.read_text(errors="ignore"))
        except OSError:
            pista = ""
        if pista:
            res.error += f" — {pista}"
        return res
    try:
        res.result = qeout.read_xml(str(job.directory))
        res.ok = res.result is not None and res.result.converged is True
        if res.result is not None and res.result.converged is False:
            # Terminó, pero sin converger. Se marca: la energía existe y es
            # un número perfectamente formado que no significa nada.
            res.ok = False
            res.error = ("finished WITHOUT CONVERGING (exhausted electron_maxstep). "
                         "Raise conv_thr, lower mixing_beta or use "
                         "'olla-dft doctor' on this folder")
        elif res.result is not None and res.result.converged is None:
            res.error = ("finished, but the XML does not confirm convergence; "
                         "the result will not be taken as valid")
    except Exception as exc:                                # noqa: BLE001
        res.error = f"finished but the result could not be read: {exc}"
    return res


def _linea_resultado(i: int, total: int, r: JobResult) -> str:
    cabeza = f"  [{i:>{len(str(total))}d}/{total}] {r.job.name} ... "
    if r.skipped:
        return cabeza + "already done"
    if r.ok:
        e = r.energy
        extra = f"  E = {e / qeout.RY_EV:.6f} Ry" if e is not None else ""
        return cabeza + f"{r.seconds:.1f} s{extra}"
    return cabeza + f"FAILED ({r.error})"


def run_all(jobs: list, pw_cmd: str = None, nproc: int = None,
            timeout: float = None, verbose: bool = True,
            paralelo: int = 1, rehacer: bool = False,
            presupuesto: float = None) -> list:
    """Corre todos los cálculos y devuelve la lista de JobResult, en orden.

    `paralelo` es cuántos cálculos simultáneos; los procesos MPI de cada uno
    salen de `reparto()` si no se fija `nproc`.

    `presupuesto` es un límite de tiempo TOTAL en segundos. Al agotarse no se
    lanzan cálculos nuevos, pero los que ya están corriendo se dejan
    terminar: matarlos a media iteración deja la carpeta en un estado del
    que no se puede reanudar, que es justo lo contrario de lo que se busca.
    """
    par, n, aviso = reparto(paralelo, nproc)
    cmd = build_command(pw_cmd, n)
    found = check_available(pw_cmd, n)
    # resolver_ejecutable también busca fuera del PATH (Homebrew, /opt/qe,
    # etc.); hay que ejecutar esa ruta, no el nombre original que fallaría.
    cmd[_qe_executable_index(cmd)] = found
    total = len(jobs)
    resultados = [None] * total
    t0 = time.time()
    candado = threading.Lock()
    agotado = threading.Event()

    if verbose:
        cab = f"Running {total} calculations with: {' '.join(cmd)}"
        if par > 1:
            cab += (f"\n  {par} at a time, {n} process{'es' if n > 1 else ''} "
                    f"each ({nucleos()} threads available). "
                    f"They finish out of order.")
        print(cab)
        if aviso:
            print(f"  WARNING: {aviso}")
        libre = memoria_libre_gb()
        if par > 1 and libre is not None and libre / par < 1.0:
            print(f"  WARNING: {libre:.1f} GB free remain, that is "
                  f"{libre / par:.1f} GB per calculation. If pw.x runs out of "
                  f"memory the system starts swapping and everything stalls; "
                  f"lower --jobs.")
        if presupuesto:
            txt = (f"{presupuesto:.0f} s" if presupuesto < 120
                   else f"{presupuesto / 60:.0f} min" if presupuesto < 7200
                   else f"{presupuesto / 3600:.1f} h")
            print(f"  Budget: {txt}. Once exhausted no more are launched, and "
                  f"those already running finish.")

    def _uno(idx_job):
        idx, job = idx_job
        if agotado.is_set():
            r = JobResult(job=job)
            r.error = "not launched: the time budget ran out"
            return idx, r
        r = run_one(job, cmd, timeout=timeout, rehacer=rehacer)
        if presupuesto and (time.time() - t0) >= presupuesto:
            agotado.set()
        if verbose:
            with candado:
                print(_linea_resultado(idx + 1, total, r), flush=True)
        return idx, r

    if par > 1 and total > 1:
        with ThreadPoolExecutor(max_workers=par) as pool:
            for idx, r in pool.map(_uno, list(enumerate(jobs))):
                resultados[idx] = r
    else:
        for idx, job in enumerate(jobs):
            _, r = _uno((idx, job))
            resultados[idx] = r

    results = [r for r in resultados if r is not None]
    if verbose:
        transcurrido = time.time() - t0
        hechos = [r for r in results if r.ok]
        sin_lanzar = [r for r in results if "time budget" in (r.error or "")]
        bad = [r for r in results if not r.ok and r not in sin_lanzar]
        tt = (f"{transcurrido:.0f} s" if transcurrido < 120
              else f"{transcurrido / 60:.1f} min" if transcurrido < 7200
              else f"{transcurrido / 3600:.2f} h")
        print(f"\nTotal time: {tt} ({len(hechos)} of {total} succeeded)")
        if sin_lanzar:
            print(f"\n{len(sin_lanzar)} calculations were NOT launched because of "
                  f"the time budget:")
            for r in sin_lanzar[:8]:
                print(f"  {r.job.name}")
            if len(sin_lanzar) > 8:
                print(f"  ... and {len(sin_lanzar) - 8} more")
            print("Run the same command again: those already finished are "
                  "skipped automatically.")
        if bad:
            print(f"\n{len(bad)} of {total} calculations failed:")
            for r in bad:
                print(f"  {r.job.name}: {r.error}")
            print("Check the corresponding outputs; the rest of the analysis "
                  "continues with the points that did finish.")
    return results


def collect(path: str, pattern: str = "*", prefix: str = None) -> list:
    """Lee los resultados ya existentes de un barrido hecho por separado.

    Sirve cuando el usuario corrió los cálculos por su cuenta (sin --run):
    recorre las subcarpetas y devuelve [(carpeta, QEResult), ...].
    """
    base = Path(path)
    found = []
    for d in sorted(base.glob(pattern)):
        if not d.is_dir():
            continue
        try:
            found.append((d, qeout.read_xml(str(d), prefix)))
        except Exception:
            continue
    return found


# ----------------------------------------------------------------------
# Diagnóstico de una corrida fallida
# ----------------------------------------------------------------------
#: Causas frecuentes de que un binario de QE ni siquiera arranque, con lo
#: que hay que hacer. Se buscan en la salida tal cual; el orden importa
#: porque la primera que empata es la que se reporta.
CAUSAS = (
    ("attempt to run as root",
     "mpirun refuses to run as root. Run as a normal user, or "
     "export OMPI_ALLOW_RUN_AS_ROOT=1 and OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1."),
    ("command not found",
     "the executable is not in the PATH. Check 'olla-dft config show' and "
     "point pw_cmd to the correct binary."),
    ("cannot open file",
     "a file it expected was not found: almost always the pseudopotential "
     "or the outdir folder. Check pseudo_dir and that the previous calculation "
     "left its .save."),
    ("reading pseudopotential",
     "the pseudopotential could not be read: misspelled name, truncated "
     "file or a format this QE does not support."),
    ("buffer overflow detected",
     "pw.x aborted inside its own binary. If it also happens with a "
     "known input, that Quantum ESPRESSO build is incompatible "
     "with the current system: install a more recent version or recompile it. "
     "If it only happens with one UPF, change that pseudopotential."),
    ("wrong ibrav",
     "the cell does not match the declared ibrav."),
    ("charge is wrong",
     "the charge does not match the valence electrons of the "
     "pseudopotentials: almost always an atom is missing or extra."),
    ("S matrix not positive definite",
     "the basis is ill-conditioned: raise ecutwfc or check for atoms "
     "that are too close."),
    ("Not enough space allocated for radial FFT",
     "raise ecutrho (8-12 times ecutwfc is usually enough for "
     "ultrasoft pseudopotentials)."),
    ("out of memory",
     "it ran out of memory. Reduce the k mesh, use fewer MPI processes or "
     "shrink the supercell."),
    ("k-point algorithm is not tested",
     "TDDFPT only implements the gamma case: the scf must use "
     "K_POINTS gamma, which is NOT the same as a 1x1x1 mesh."),
    ("Linear response calculation" ,
     "TDDFPT does not support symmetry: the previous scf must use "
     "nosym=.true. and noinv=.true."),
    ("some of the original symmetry operations not satisfied",
     "the atoms moved and broke the symmetry that pw.x detected at the "
     "start. In molecular dynamics and in relaxations from a symmetric "
     "configuration you must set nosym=.true."),
    ("too many bands are not converged",
     "the diagonalizer does not converge: lower mixing_beta, raise "
     "electron_maxstep, or try diagonalization='cg'."),
    ("SCF correction compared to forces is large",
     "the scf is not converged tightly enough for the forces being "
     "requested: lower conv_thr."),
)


def failure_hint(text: str) -> str:
    """Traduce la salida de un binario de QE a una causa probable, si la hay."""
    for clave, explicacion in CAUSAS:
        if clave.lower() in text.lower():
            return explicacion
    return ""


def failure_message(stem: str, out_file, text: str = None,
                    lineas: int = 6) -> str:
    """Mensaje de fallo con el final del log y, si se reconoce, la causa.

    Decir solo "revisa tal archivo" obliga a abrirlo para enterarse de que
    faltaba un pseudopotencial. Las últimas líneas del log casi siempre
    traen el motivo, así que se muestran aquí.
    """
    out_file = Path(out_file)
    if text is None:
        try:
            text = out_file.read_text(errors="ignore")
        except OSError:
            text = ""
    partes = [f"{stem} failed; check {out_file}"]
    pista = failure_hint(text)
    if pista:
        partes.append(f"Probable cause: {pista}")
    cola = [ln.rstrip() for ln in text.splitlines() if ln.strip()]
    # Si QE alcanzó a escribir su propio bloque de error, ese es el que
    # importa; si no, las últimas líneas de lo que haya.
    marca = next((i for i, ln in enumerate(cola)
                  if "Error in routine" in ln), None)
    if marca is not None:
        cola = cola[marca:marca + lineas]
    else:
        cola = cola[-lineas:]
    if cola:
        partes.append("End of the log:")
        partes += [f"  {ln}" for ln in cola]
    return "\n".join(partes)
