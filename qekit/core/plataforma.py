# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Dónde van las cosas y cómo se llaman en cada sistema operativo.

Olla-DFT se escribió y se probó en Linux. Casi todo es portable por accidente
—Python, numpy, ASE lo son—, pero hay cuatro cosas que no lo son y que hay
que decidir a propósito:

  1. **Dónde vive la configuración.** `~/.config/qekit` es la convención de
     Linux (XDG). En macOS lo suyo es `~/Library/Application Support`, y en
     Windows `%APPDATA%`. Escribir en `~/.config` en Windows funciona, pero
     deja una carpeta oculta donde nadie la busca y que ningún desinstalador
     limpia.
  2. **Cómo se llaman los binarios.** En Windows los ejecutables de Quantum
     ESPRESSO son `pw.exe`, no `pw.x`. Buscar solo `pw.x` hace que Olla-DFT diga
     "no encuentro pw.x" en una máquina donde está instalado.
  3. **Cómo se lanza MPI.** `mpirun` en OpenMPI, `mpiexec` en MPICH y en MS-MPI,
     `srun` en un clúster con Slurm. El valor por defecto de Olla-DFT era
     `mpirun`, que en Windows no existe.
  4. **Cómo se cuentan los núcleos en un guion de shell.** `nproc` es de
     coreutils y en macOS no está: allí es `sysctl -n hw.ncpu`.

Este módulo concentra esas cuatro decisiones para que el resto del código no
tenga que preguntar en qué sistema está.

Lo que NO hace este módulo es fingir que Quantum ESPRESSO es portable. En
Windows hay tres caminos reales —WSL, binarios nativos y Cygwin— y no dan la
misma experiencia; `olla-dft sistema` los explica con lo que encuentre en la
máquina concreta.
"""

import os
import shutil
import sys
from pathlib import Path

WINDOWS = sys.platform.startswith("win")
MACOS = sys.platform == "darwin"
LINUX = sys.platform.startswith("linux")

#: Nombre del sistema, para informes
def nombre_sistema() -> str:
    if WINDOWS:
        return "Windows"
    if MACOS:
        return "macOS"
    if LINUX:
        return "Linux"
    return sys.platform


APP_DIR = "olla-dft"        # nombre de la carpeta de configuración y datos
_APP_DIR_HEREDADO = "QEkit"  # el que usaban las versiones anteriores a la 1.0


def _env_dir(*nombres) -> Path:
    for nombre in nombres:
        forzado = os.environ.get(nombre)
        if forzado:
            return Path(forzado).expanduser()
    return None


def dir_config() -> Path:
    """Carpeta de configuración, según la convención de cada sistema.

    Se respeta ``OLLA_DFT_CONFIG_DIR`` (o el heredado ``QEKIT_CONFIG_DIR``)
    por encima de todo: es lo que permite correr Olla-DFT desde un lápiz USB
    o en un clúster con el HOME lleno, y lo que hace que las pruebas no
    toquen la configuración de verdad.
    """
    forzado = _env_dir("OLLA_DFT_CONFIG_DIR", "QEKIT_CONFIG_DIR")
    if forzado:
        return forzado
    if WINDOWS:
        base = os.environ.get("APPDATA")
        if base:
            return Path(base) / APP_DIR
        return Path.home() / "AppData" / "Roaming" / APP_DIR
    if MACOS:
        return Path.home() / "Library" / "Application Support" / APP_DIR
    xdg = os.environ.get("XDG_CONFIG_HOME")
    return (Path(xdg) if xdg else Path.home() / ".config") / APP_DIR


def dir_data() -> Path:
    """Carpeta para datos grandes y reemplazables, según cada plataforma.

    Los modelos y el histórico no son configuración: guardarlos en
    ``~/.config`` llena copias de seguridad y rompe las convenciones de cada
    escritorio. ``OLLA_DFT_DATA_DIR`` (o ``QEKIT_DATA_DIR``) permite, además,
    ponerlos en un disco con más espacio.
    """
    forced = _env_dir("OLLA_DFT_DATA_DIR", "QEKIT_DATA_DIR")
    if forced:
        return forced
    if WINDOWS:
        base = os.environ.get("LOCALAPPDATA") or os.environ.get("APPDATA")
        if base:
            return Path(base) / APP_DIR
        return Path.home() / "AppData" / "Local" / APP_DIR
    if MACOS:
        return Path.home() / "Library" / "Application Support" / APP_DIR
    xdg = os.environ.get("XDG_DATA_HOME")
    return (Path(xdg) if xdg else Path.home() / ".local" / "share") / APP_DIR


def dirs_config_heredados() -> list:
    """Donde guardaban la configuración las versiones anteriores, para migrarla.

    Hasta la 0.27 siempre ``~/.config/qekit``; de la 0.28 a la 0.35, la
    carpeta ``QEkit`` de cada sistema.
    """
    candidatos = [Path.home() / ".config" / "qekit"]
    if WINDOWS:
        base = os.environ.get("APPDATA")
        candidatos.append((Path(base) if base else Path.home() / "AppData"
                           / "Roaming") / _APP_DIR_HEREDADO)
    elif MACOS:
        candidatos.append(Path.home() / "Library" / "Application Support"
                          / _APP_DIR_HEREDADO)
    return candidatos


def dir_config_heredada() -> Path:
    """La carpeta heredada más antigua (compatibilidad)."""
    return dirs_config_heredados()[0]


#: Los ejecutables de Quantum ESPRESSO, con las dos terminaciones posibles.
def nombres_ejecutable(base: str) -> list:
    """'pw' -> ['pw.x', 'pw.exe'] en Windows, ['pw.x'] en el resto.

    El orden importa: se prueba primero el nombre nativo del sistema.
    """
    base = str(base)
    if base.endswith((".x", ".exe")):
        base = base.rsplit(".", 1)[0]
    if WINDOWS:
        return [f"{base}.exe", f"{base}.x"]
    return [f"{base}.x", f"{base}.exe"]


def buscar_ejecutable(base: str, extra_dirs=()) -> str:
    """Devuelve la ruta del binario si está, o None.

    Busca por PATH y, además, en las carpetas que se le pasen: en macOS con
    Homebrew y en muchas instalaciones a mano, los binarios de QE no están en
    el PATH.
    """
    for nombre in nombres_ejecutable(base):
        hallado = shutil.which(nombre)
        if hallado:
            return hallado
        for d in extra_dirs:
            cand = Path(d).expanduser() / nombre
            if cand.exists():
                return str(cand)
    return None


#: Carpetas donde suele estar QE en cada sistema, para buscar sin PATH.
def dirs_probables_qe() -> list:
    if WINDOWS:
        pf = os.environ.get("ProgramFiles", r"C:\Program Files")
        return [Path(pf) / "QE" / "bin", Path(pf) / "quantum-espresso" / "bin",
                Path("C:/qe/bin")]
    if MACOS:
        return [Path("/opt/homebrew/bin"), Path("/usr/local/bin"),
                Path("/opt/local/bin"),              # MacPorts
                Path.home() / "q-e" / "bin"]
    return [Path("/usr/bin"), Path("/usr/local/bin"),
            Path("/opt/qe/bin"), Path.home() / "q-e" / "bin",
            Path("/usr/lib64/openmpi/bin")]


def lanzador_mpi() -> str:
    """La plantilla de lanzamiento MPI que exista en esta máquina.

    Se prueban en orden de probabilidad. Si no hay ninguno, se devuelve
    cadena vacía: Olla-DFT corre en serie, que es lo correcto en un portátil.
    """
    for cand, plantilla in (("mpirun", "mpirun -np {n}"),
                            ("mpiexec", "mpiexec -n {n}"),
                            ("srun", "srun -n {n}")):
        if shutil.which(cand):
            return plantilla
    return ""


def cuenta_nucleos_shell() -> str:
    """Trozo de shell POSIX que cuenta núcleos en Linux Y en macOS.

    `nproc` es de coreutils y en macOS no está; allí es `sysctl -n hw.ncpu`.
    Un guion generado con solo `nproc` deja a los usuarios de Mac con un
    número vacío y una división por cero.
    """
    return "$( (nproc 2>/dev/null) || (sysctl -n hw.ncpu 2>/dev/null) || echo 1 )"


def escribir_script(ruta, texto: str, ejecutable: bool = True) -> Path:
    """Escribe un guion con finales de línea POSIX y permiso de ejecución.

    Los dos detalles importan. Si Python abre el archivo en modo texto en
    Windows, escribe CRLF, y un `.sh` con CRLF ejecutado en WSL o en un
    clúster falla con «bad interpreter: /bin/bash^M», que es un error
    célebre por lo difícil de leer que es. Y `chmod` no existe en Windows:
    hay que intentarlo sin que reviente.
    """
    p = Path(ruta)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", encoding="utf-8", newline="\n") as f:
        f.write(texto)
    if ejecutable and not WINDOWS:
        try:
            p.chmod(0o755)
        except OSError:
            pass
    return p


_SECUENCIAL_PY = r'''#!/usr/bin/env python3
# Runs the steps of this calculation, in order. Generated by Olla-DFT.
# Does the same as correr.sh and works equally on Windows, where there is no bash.
import os
import shutil
import subprocess
import sys
from pathlib import Path

PASOS = __PASOS__          # (executable, input, output)
NPROC_DEF = __NPROC__
AQUI = Path(__file__).resolve().parent


def binario(nombre):
    base = nombre
    for suf in (".x", ".exe"):
        if base.endswith(suf):
            base = base[: -len(suf)]
    for cand in (nombre, base + ".exe", base + ".x"):
        hallado = shutil.which(cand)
        if hallado:
            return hallado
        if Path(cand).exists():
            return str(Path(cand).resolve())
    return None


def main():
    np_ = int(os.environ.get("NPROC", str(NPROC_DEF)))
    lanz = []
    if np_ > 1:
        for cand, bandera in (("mpirun", "-np"), ("mpiexec", "-n")):
            if shutil.which(cand):
                lanz = [cand, bandera, str(np_)]
                break
    for exe, entrada, salida in PASOS:
        ruta = binario(exe)
        if ruta is None:
            print("Cannot find '%s'. Install it or put it in the PATH." % exe)
            return 2
        print(">> %s < %s" % (exe, entrada), flush=True)
        with open(AQUI / entrada) as fi, open(AQUI / salida, "w") as fo:
            r = subprocess.run(lanz + [ruta], stdin=fi, stdout=fo,
                               stderr=subprocess.STDOUT, cwd=str(AQUI))
        if r.returncode != 0:
            print("  failure in %s (code %d). See %s." %
                  (exe, r.returncode, salida))
            return 1
    print("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
'''


def build_sequential_python_script(pasos, nproc: int = 1) -> str:
    """Construye el guion Python portable para una cadena de pasos de QE."""
    return (_SECUENCIAL_PY
            .replace("__PASOS__", repr([tuple(p) for p in pasos]))
            .replace("__NPROC__", str(max(1, int(nproc)))))


def escribir_par_de_guiones(carpeta, pasos, nombre="correr") -> list:
    """Escribe `<nombre>.sh` y `<nombre>.py` para una cadena secuencial.

    `pasos` es una lista de (ejecutable, entrada, salida). Los dos guiones
    hacen lo mismo; el `.py` es el que sirve en Windows y el que sabe que el
    binario puede llamarse `.exe`.

    Existe porque media docena de módulos escribían su propio `correr.sh` a
    mano, cada uno con su `#!/bin/bash` y su `pw.x` cableado. Con eso, en
    Windows no había forma de correr una cadena de fonones o de Wannier
    aunque Quantum ESPRESSO estuviera instalado.
    """
    from pathlib import Path as _P
    carpeta = _P(carpeta)
    sh = ["#!/bin/bash", "# Generated by Olla-DFT — the steps, in order.",
          "# On Windows, or without bash:  python %s.py" % nombre, "set -e",
          'NP=${NPROC:-1}',
          'if [ "$NP" -gt 1 ] && command -v mpirun >/dev/null 2>&1; then',
          '  LANZ="mpirun -np $NP"',
          'elif [ "$NP" -gt 1 ] && command -v mpiexec >/dev/null 2>&1; then',
          '  LANZ="mpiexec -n $NP"', "else", '  LANZ=""', "fi", ""]
    for exe, entrada, salida in pasos:
        sh.append('echo ">> %s < %s"' % (exe, entrada))
        sh.append('$LANZ %s -in %s > %s 2>&1' % (exe, entrada, salida))
    sh.append('echo "Done."')
    escribir_script(carpeta / f"{nombre}.sh", "\n".join(sh) + "\n")
    escribir_script(carpeta / f"{nombre}.py",
                    build_sequential_python_script(pasos))
    return [str(carpeta / f"{nombre}.sh"), str(carpeta / f"{nombre}.py")]


def informe() -> str:
    """Qué ve Olla-DFT de esta máquina, y qué hacer si algo falta.

    Es el primer comando que hay que correr en un sistema nuevo. Contesta a
    la vez las cuatro preguntas que llevan a abandonar un programa
    científico: dónde guarda mis cosas, encuentra mi Quantum ESPRESSO,
    puede escribir los caracteres que usa, y cómo lanzo los cálculos aquí.
    """
    import platform
    import sys as _sys
    from qekit import __version__
    from qekit.core import consola

    L = [f"--- Olla-DFT {__version__} on this machine ---",
         f"System:    {nombre_sistema()}  ({platform.platform()})",
         f"Python:    {_sys.version.split()[0]}  ({_sys.executable})",
         ""]

    # --- salida ---
    est = consola.ESTADO
    cod = est.get("codificacion") or getattr(_sys.stdout, "encoding", "?")
    L += ["Text output", f"  encoding: {cod}"]
    if est.get("forzado_utf8"):
        L.append("  UTF-8 was forced on a console that did not have it enabled.")
    if est.get("translitera"):
        L += ["  this console does NOT support Å, α or →: they are transliterated to ASCII.",
              "  On Windows 10 or later you can enable UTF-8 like this:",
              "      chcp 65001",
              "  or set PYTHONUTF8=1 in the environment variables."]
    else:
        L.append("  supports Å α → ① ✓ without problems.")
    L.append("  (with --ascii the output is forced to ASCII in any case)")

    # --- configuración ---
    d = dir_config()
    L += ["", "Configuration",
          f"  folder:    {d}",
          f"  file:      {'exists' if (d / 'config.ini').exists() else 'not yet (created on the first `olla-dft config set`)'}"]
    vieja = dir_config_heredada()
    if vieja != d and (vieja / "config.ini").exists():
        L.append(f"  there is an old configuration in {vieja}; it is copied "
                 f"automatically the first time.")
    L.append("  it can be moved with the OLLA_DFT_CONFIG_DIR variable.")

    # --- Quantum ESPRESSO ---
    L += ["", "Quantum ESPRESSO"]
    dirs = dirs_probables_qe()
    hallados, faltan = [], []
    for base in ("pw", "ph", "dos", "projwfc", "bands", "pp", "epsilon",
                 "q2r", "matdyn", "pw2wannier90", "hp", "neb", "ld1",
                 "xspectra", "turbo_lanczos", "pwcond"):
        ruta = buscar_ejecutable(base, dirs)
        (hallados if ruta else faltan).append((base, ruta))
    if hallados:
        L.append(f"  found {len(hallados)} of "
                 f"{len(hallados) + len(faltan)} executables.")
        L.append(f"  pw:  {hallados[0][1]}")
    else:
        L += ["  NONE were found.",
              "  Olla-DFT is still useful for generating inputs and for "
              "post-processing outputs",
              "  you bring from another machine; only --run needs the "
              "binaries here."]
    if faltan and hallados:
        L.append("  missing: " + ", ".join(b for b, _ in faltan[:8])
                 + ("..." if len(faltan) > 8 else ""))
        L.append("  (that is fine: each module warns if it needs one that is "
                 "not there)")

    # --- paralelismo ---
    L += ["", "Parallelism"]
    lanz = lanzador_mpi()
    L.append(f"  visible cores: {_nucleos()}")
    if lanz:
        L.append(f"  MPI launcher:  {lanz.format(n=_nucleos())}")
    else:
        L += ["  no mpirun/mpiexec: calculations will run in serial.",
              "  On a laptop that is usually the right thing; running several "
              "points at once",
              "  with `-j N` makes better use of 2 or 4 cores than a single pw.x "
              "with MPI."]

    L += ["", "How to launch calculations here"]
    L += ["  " + x for x in _consejo_plataforma()]
    return "\n".join(L)


def _nucleos() -> int:
    import os as _os
    try:
        return len(_os.sched_getaffinity(0))
    except AttributeError:
        return _os.cpu_count() or 1


def _consejo_plataforma() -> list:
    if WINDOWS:
        return [
            "Olla-DFT works on Windows: it is pure Python. Quantum ESPRESSO is",
            "another matter, and there are three routes, from most to least recommended:",
            "",
            "  1. WSL2 (the simplest and the best tested).",
            "     Install Ubuntu with `wsl --install`, and inside it:",
            "         sudo apt install quantum-espresso python3-pip",
            "         pip install olla_dft-*.whl",
            "     You work inside WSL and everything behaves as on Linux.",
            "",
            "  2. Native QE binaries for Windows (they are called pw.exe).",
            "     Olla-DFT looks for them itself; if they are elsewhere:",
            "         olla-dft config set pw_cmd C:\\\\path\\\\to\\\\pw.exe",
            "     The generated scripts are run with `python run.py`, not",
            "     with run.sh: there is no bash on Windows.",
            "",
            "  3. Olla-DFT here and Quantum ESPRESSO on a cluster.",
            "     Generate the inputs on your laptop, copy them, run them there and",
            "     bring back the outputs: all post-processing works without",
            "     a single QE binary on this machine.",
            "",
            "And if the console shows strange symbols:  chcp 65001",
        ]
    if MACOS:
        return [
            "Quantum ESPRESSO is installed with Homebrew or with MacPorts:",
            "     brew install quantum-espresso",
            "     sudo port install quantum-espresso",
            "On Apple Silicon, Homebrew installs into /opt/homebrew/bin, which is not",
            "always in the PATH of a non-interactive shell. Olla-DFT looks",
            "there anyway.",
            "",
            "The generated scripts use `sysctl -n hw.ncpu` to count",
            "cores, not `nproc`, which does not exist on macOS.",
        ]
    return [
        "Quantum ESPRESSO is packaged in the major distributions:",
        "     Debian/Ubuntu:  sudo apt install quantum-espresso",
        "     Fedora/RHEL:    sudo dnf install quantum-espresso",
        "     openSUSE:       sudo zypper install quantum-espresso",
        "     Arch:           sudo pacman -S quantum-espresso",
        "     conda:          conda install -c conda-forge qe",
        "If you compiled it by hand and it is not in the PATH:",
        "     olla-dft config set pw_cmd /path/to/bin/pw.x",
    ]
