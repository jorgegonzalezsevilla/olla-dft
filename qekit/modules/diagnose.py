# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Diagnóstico de un cálculo de pw.x: ¿sirve, y si no, por qué?

Todo lo que este módulo usa ya estaba en los archivos que QE deja: el XML
guarda si convergió, en cuántos pasos y con qué error; el stdout guarda la
historia completa de iteraciones SCF y, en un relax, cada paso iónico.
Nadie los mira, y son justo los que dicen si el resultado sirve.

LO QUE DISTINGUE ESTE MÓDULO
----------------------------
Un SCF que no converge tiene al menos dos causas con remedios OPUESTOS:

- **oscilación de carga** (charge sloshing): el error sube y baja en vez de
  bajar. Típico de losas grandes, metales y celdas con vacío. El remedio es
  mezclar MENOS (bajar mixing_beta) y usar mixing_mode='local-TF';
- **convergencia lenta monótona**: el error baja siempre, pero demasiado
  despacio. Ahí el remedio es el contrario, mezclar MÁS (subir
  mixing_beta) o simplemente dar más pasos.

Aplicar el remedio equivocado empeora el problema, así que el módulo
distingue los dos casos por la forma de la curva en vez de dar un consejo
genérico.
"""

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from qekit.core import qeout
from qekit.core import style as qstyle

_RE_ITER = re.compile(r"iteration #\s*(\d+)\s+ecut=\s*([\d.]+)\s*Ry\s+beta=\s*([\d.]+)")
_RE_ACC = re.compile(r"estimated scf accuracy\s*<\s*([\dEe.+-]+)\s*Ry")
_RE_ETOT = re.compile(r"total energy\s*=\s*(-?[\dEe.+-]+)\s*Ry")
_RE_FORCE = re.compile(r"Total force\s*=\s*([\dEe.+-]+)")
_RE_PRESS = re.compile(r"P=\s*(-?[\dEe.+-]+)")
_RE_WARN = re.compile(r"^\s*(Warning|WARNING|%%%%)", re.M)
_RE_MAXSTEP = re.compile(r"convergence NOT achieved after\s*(\d+)\s*iterations")


@dataclass
class ScfHistory:
    accuracy: list = field(default_factory=list)     # Ry, por iteración (último ciclo)
    energies: list = field(default_factory=list)     # Ry (último ciclo)
    beta: float = None
    converged: bool = None
    n_iter: int = 0
    n_ciclos: int = 1            # ciclos SCF vistos (1 en un scf, N en un relax)
    patologia: str = ""          # "oscilacion" | "lenta" | "estancada" | ""
    consejo: str = ""


@dataclass
class Trajectory:
    energies: list = field(default_factory=list)     # Ry por paso iónico
    forces: list = field(default_factory=list)       # Ry/bohr (Total force)
    pressures: list = field(default_factory=list)    # kbar
    n_steps: int = 0


@dataclass
class Diagnosis:
    result: qeout.QEResult = None
    scf: ScfHistory = None
    traj: Trajectory = None
    warnings: list = field(default_factory=list)
    problemas: list = field(default_factory=list)
    stdout_path: str = ""


def _ciclos_scf(texto: str) -> list:
    """Trozos del stdout, uno por ciclo SCF.

    En un relax cada paso iónico arranca su propio ciclo con
    'iteration #  1'; se corta ahí. Un scf da un solo trozo.
    """
    cortes = [m.start() for m in _RE_ITER.finditer(texto)
              if int(m.group(1)) == 1]
    if not cortes:
        return [texto]
    cortes.append(len(texto))
    return [texto[a:b] for a, b in zip(cortes, cortes[1:])]


def read_scf_history(stdout_path) -> ScfHistory:
    """Historia de iteraciones SCF del stdout de pw.x, y su diagnóstico.

    En un relax hay un ciclo SCF por paso iónico y son independientes:
    concatenar sus errores mezclaría el final convergido de uno con el
    arranque del siguiente y parecería una oscilación. Se clasifica solo
    el ÚLTIMO ciclo, que es el que decide si el cálculo terminó bien, y
    `n_ciclos` dice cuántos se vieron.
    """
    texto = Path(stdout_path).read_text(errors="ignore")
    h = ScfHistory()
    ciclos = _ciclos_scf(texto)
    ultimo = ciclos[-1]
    h.n_ciclos = len(ciclos)
    h.accuracy = [float(x) for x in _RE_ACC.findall(ultimo)]
    h.energies = [float(x) for x in _RE_ETOT.findall(ultimo)]
    betas = _RE_ITER.findall(texto)
    if betas:
        h.beta = float(betas[0][2])
    h.n_iter = len(h.accuracy)
    h.converged = "convergence has been achieved" in ultimo
    if _RE_MAXSTEP.search(ultimo):
        h.converged = False
    _clasificar(h)
    return h


def _clasificar(h: ScfHistory) -> None:
    """Distingue oscilacion de convergencia lenta por la FORMA de la curva.

    Se miran DOS cosas, no una: con que frecuencia sube el error y CUANTO
    sube. Unas pocas subidas enormes son mejor senal de oscilacion que
    muchas subidas minimas, asi que contar solo la frecuencia se equivoca.

    Las dos primeras iteraciones se ignoran: un salto grande al arranque es
    un transitorio normal mientras la densidad inicial se acomoda, no una
    patologia.
    """
    a = np.array(h.accuracy, dtype=float)
    if a.size < 2 or h.converged:
        return

    beta = h.beta if h.beta is not None else 0.4

    # Con muy pocas iteraciones no hay forma honesta de distinguir
    # oscilacion de lentitud: una sola subida entre tres diferencias ya da
    # 33 %, y eso es ruido, no diagnostico. Mejor decirlo que inventarlo.
    if a.size < 8:
        h.patologia = "pocos_datos"
        h.consejo = (
            f"only {a.size} iterations: not enough to distinguish "
            "charge sloshing from\nslow convergence, which call for opposite "
            "remedies. Raise electron_maxstep\n(to 100 or more) and look at "
            "the curve again; with the cycle cut so early,\nany "
            "diagnosis would be guesswork.")
        return

    cola = a[2:] if a.size > 5 else a          # saltar el transitorio
    difs = np.diff(cola)
    subidas = int(np.sum(difs > 0))
    frac_subidas = subidas / max(len(cola) - 1, 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        razones = cola[1:] / np.maximum(cola[:-1], 1e-300)
        peor = float(np.max(razones)) if razones.size else 1.0
        decadas = np.log10(a[0] / a[-1]) if a[-1] > 0 else np.inf

    # el criterio de frecuencia solo se aplica con suficientes puntos
    frecuente = len(cola) >= 6 and frac_subidas > 0.25
    if frecuente or peor > 5.0:
        motivo = (f"the error rose in {subidas} of {len(cola)-1} "
                  "iterations" if frecuente
                  else f"the error grew by a factor of {peor:.0f} in one "
                       "iteration")
        h.patologia = "oscilacion"
        h.consejo = (
            f"{motivo}: this is CHARGE SLOSHING, typical of slabs,\n"
            "metals and cells with a lot of vacuum. It is fixed by mixing LESS:\n"
            f"  mixing_beta = {max(0.05, beta / 3):.2f}   "
            f"(now {beta:.2f})\n"
            "  mixing_mode = 'local-TF'   (designed precisely for this case)\n"
            "  mixing_ndim = 12           (more mixing history)\n"
            "Raising mixing_beta here makes it WORSE.")
    elif decadas < 3:
        h.patologia = "estancada"
        h.consejo = (
            f"in {len(a)} iterations the error only dropped {decadas:.1f} "
            "orders of magnitude and\nwent flat: this is not slowness, it is "
            "stalled. It is usually a badly posed magnetic\nor occupation "
            "state, or a structure with atoms almost on top of each other.\nCheck "
            "starting_magnetization, the smearing and the interatomic "
            "distances.")
    else:
        h.patologia = "lenta"
        # Si beta ya es agresivo, subirlo mas seria justo lo contrario de
        # lo que conviene: ahi lo que falta son pasos, no mezcla.
        if beta >= 0.6:
            h.consejo = (
                f"the error decreases monotonically and reached {a[-1]:.1e} Ry "
                "without reaching the\nthreshold: it lacked STEPS, not mixing. "
                f"mixing_beta is already at {beta:.2f}, which is\naggressive; "
                "raising it further risks destabilizing it.\n"
                "  electron_maxstep = 300")
        else:
            h.consejo = (
                "the error decreases monotonically but did not reach the threshold: "
                "this is SLOW\nconvergence, not sloshing. Here it does help to "
                "mix MORE or allow more steps:\n"
                f"  mixing_beta = {min(0.7, max(beta * 1.75, 0.3)):.2f}   "
                f"(now {beta:.2f})\n"
                "  electron_maxstep = 300")


def read_trajectory(stdout_path) -> Trajectory:
    """Pasos iónicos de un relax/vc-relax desde el stdout."""
    texto = Path(stdout_path).read_text(errors="ignore")
    t = Trajectory()
    # las energías de cada paso convergido llevan '!' delante
    t.energies = [float(m) for m in re.findall(
        r"^!\s+total energy\s*=\s*(-?[\dEe.+-]+)\s*Ry", texto, re.M)]
    t.forces = [float(x) for x in _RE_FORCE.findall(texto)]
    t.pressures = [float(x) for x in _RE_PRESS.findall(texto)]
    t.n_steps = len(t.energies)
    return t


def find_stdout(workdir) -> Path:
    """Busca el stdout de pw.x en una carpeta de cálculo."""
    d = Path(workdir)
    if d.is_file():
        return d
    for patron in ("*.out", "pw.out", "scf.out", "relax.out", "out.*"):
        for f in sorted(d.glob(patron)):
            try:
                head = f.read_text(errors="ignore")[:4000]
            except OSError:
                continue
            if "Program PWSCF" in head:
                return f
    return None


def diagnose(workdir, prefix: str = None) -> Diagnosis:
    """Lee XML y stdout de una carpeta y arma el diagnóstico."""
    d = Diagnosis()
    try:
        xml = qeout.find_xml(str(workdir), prefix)
        d.result = qeout.read_xml(xml)
    except (FileNotFoundError, ValueError) as exc:
        d.problemas.append(f"could not read the XML: {exc}")

    so = find_stdout(workdir)
    if so is not None:
        d.stdout_path = str(so)
        d.scf = read_scf_history(so)
        d.traj = read_trajectory(so)
        texto = so.read_text(errors="ignore")
        for linea in texto.splitlines():
            if _RE_WARN.match(linea) and linea.strip() not in d.warnings:
                d.warnings.append(linea.strip())
        if "Error in routine" in texto:
            i = texto.index("Error in routine")
            d.problemas.append(texto[i:i + 200].split("%%%")[0].strip())

    r = d.result
    if r is not None:
        if r.converged is False:
            d.problemas.append("the SCF did NOT converge: the result is unusable")
        if r.max_force is not None and r.max_force > 0.05:
            d.problemas.append(
                f"maximum residual force {r.max_force:.4f} eV/Å: the "
                "structure is not relaxed (usual threshold 0.01–0.03)")
        if r.pressure is not None and abs(r.pressure) > 1.0 and \
                r.calculation in ("vc-relax", "relax", "scf"):
            d.problemas.append(
                f"residual pressure {r.pressure:+.2f} GPa: the cell is not "
                "in equilibrium with these cutoffs")
    return d


#: English display names of the `patologia` identifiers (the identifiers
#: themselves are kept as they are stored in ScfHistory).
_PATOLOGIA_EN = {"oscilacion": "charge sloshing", "lenta": "slow",
                 "estancada": "stalled", "pocos_datos": "too few iterations"}


def report(d: Diagnosis) -> str:
    r = d.result
    lines = ["--- Calculation diagnosis ---"]
    if r is not None:
        lines += [f"File: {r.xml_path}",
                  f"Type: {r.calculation or '?'}  |  "
                  f"{r.functional or '?'}  |  ecut {r.ecutwfc or '?'}/"
                  f"{r.ecutrho or '?'} Ry"]
        if r.kgrid:
            lines.append(f"k-grid: {r.kgrid[0]}x{r.kgrid[1]}x{r.kgrid[2]}"
                         f"  |  {r.nk} points in the IBZ  |  "
                         f"{r.n_sym or '?'} symmetry operations")
        estado = ("converged" if r.converged else "did NOT converge"
                  if r.converged is not None else "convergence unknown")
        extra = ""
        if r.n_scf_steps:
            extra = f" in {r.n_scf_steps} steps"
        if r.scf_error is not None:
            extra += f", final error {r.scf_error:.2e} Ry"
        lines.append(f"SCF: {estado}{extra}")
        if r.max_force is not None:
            lines.append(f"Maximum residual force: {r.max_force:.5f} eV/Å")
        if r.pressure is not None:
            lines.append(f"Residual pressure: {r.pressure:+.3f} GPa")
        if r.total_magnetization is not None and \
                abs(r.total_magnetization) > 1e-8:
            lines.append(f"Magnetization: {r.total_magnetization:.3f} μB "
                         f"(absolute {r.absolute_magnetization:.3f})")
        if r.wall_time:
            lines.append(f"Time: {r.wall_time:.1f} s wall "
                         f"({r.cpu_time:.1f} s CPU)")

    if d.traj and d.traj.n_steps > 1:
        t = d.traj
        lines += ["", f"Relaxation: {t.n_steps} ionic steps",
                  f"  energy: {t.energies[0]:.6f} -> {t.energies[-1]:.6f} Ry "
                  f"({(t.energies[-1]-t.energies[0])*13.6057:.4f} eV)"]
        if t.forces:
            lines.append(f"  total force: {t.forces[0]:.5f} -> "
                         f"{t.forces[-1]:.5f} Ry/bohr")
        subidas = sum(1 for a, b in zip(t.energies, t.energies[1:]) if b > a)
        if subidas > t.n_steps // 3:
            lines.append(
                f"  WARNING: the energy rose in {subidas} of {t.n_steps-1} "
                "steps. A healthy relaxation\n  almost always goes down; this "
                "suggests a very flat energy surface or a\n  BFGS step "
                "that is too large.")

    if d.scf and d.scf.n_iter:
        ciclos = ""
        if d.scf.n_ciclos > 1:
            ciclos = (f" in the last of {d.scf.n_ciclos} SCF cycles (one "
                      "per ionic step; only the last one is diagnosed)")
        lines += ["", f"SCF history: {d.scf.n_iter} iterations{ciclos}"
                       f"{'' if d.scf.beta is None else f', beta = {d.scf.beta:.2f}'}"]
        if d.scf.patologia:
            nombre_pat = _PATOLOGIA_EN.get(d.scf.patologia, d.scf.patologia)
            lines += ["", f"CONVERGENCE PROBLEM ({nombre_pat}):",
                      d.scf.consejo]

    if d.problemas:
        lines += ["", "PROBLEMS:"]
        lines += [f"  - {p}" for p in d.problemas]
    elif r is not None and r.converged:
        lines += ["", "No problems detected."]

    if d.warnings:
        lines += ["", f"QE warnings ({len(d.warnings)}):"]
        lines += [f"  {w}" for w in d.warnings[:6]]
    return "\n".join(lines)


def plot(d: Diagnosis, outfile: str = "diagnostico", formats="pdf,png",
         theme: str = None, family: str = None, background: str = None,
         palette=None, usetex: bool = None, width="double",
         journal: str = "generic", aspect: float = 0.42,
         mono: bool = False, dpi: int = None) -> list:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    tiene_traj = bool(d.traj and d.traj.n_steps > 1)
    n = 2 if tiene_traj else 1
    st = qstyle.apply(theme, family=family, background=background,
                      palette=palette, usetex=usetex, mono=mono)
    fig = plt.figure(figsize=qstyle.figure_size(width, journal, aspect),
                     layout="constrained")
    c = qstyle.palette(3, mono=mono)

    ax = qstyle.finish_axes(fig.add_subplot(1, n, 1))
    if d.scf and d.scf.accuracy:
        ax.semilogy(range(1, len(d.scf.accuracy) + 1), d.scf.accuracy,
                    "o-", color=c[0], lw=st["line"], ms=3)
    ax.set_xlabel("SCF iteration")
    ax.set_ylabel("estimated accuracy (Ry)")
    qstyle.panel_label(ax, "(a)")

    if tiene_traj:
        ax2 = qstyle.finish_axes(fig.add_subplot(1, n, 2))
        e = np.array(d.traj.energies) * 13.605693
        ax2.plot(range(1, len(e) + 1), e - e[-1], "o-", color=c[1],
                 lw=st["line"], ms=3)
        ax2.set_xlabel("ionic step")
        ax2.set_ylabel(r"$E - E_\mathrm{final}$ (eV)")
        qstyle.panel_label(ax2, "(b)")

    written = qstyle.save(fig, outfile, formats, dpi=dpi,
                          modulo="diagnóstico")
    plt.close(fig)
    return written
