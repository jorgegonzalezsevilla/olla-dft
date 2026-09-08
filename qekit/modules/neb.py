# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Caminos de reacción y barreras de activación con `neb.x`.

QUÉ ES
------
NEB (nudged elastic band) encuentra el camino de mínima energía entre dos
estructuras conocidas —reactivo y producto—. Interpola una cadena de
imágenes entre ellas y las relaja perpendicularmente al camino, con
muelles que las mantienen repartidas. El máximo de la curva es el estado
de transición, y su altura sobre el reactivo es la barrera.

LO QUE DECIDE SI EL RESULTADO SIRVE
-----------------------------------
1. **La imagen trepadora (CI).** Sin ella, la imagen más alta se queda
   cerca del máximo pero no encima, y la barrera sale SUBESTIMADA —
   típicamente por décimas de eV. Olla-DFT la activa por omisión y avisa si
   se desactiva.

2. **El número de imágenes.** Pocas y el máximo se pierde entre dos
   imágenes; muchas y cada paso cuesta. De 7 a 11 es lo habitual. La
   pista de que faltan: la curva interpolada tiene el máximo lejos de
   cualquier imagen calculada.

3. **Los extremos.** Reactivo y producto tienen que estar RELAJADOS con
   los mismos parámetros que se van a usar en el NEB. Si no, la barrera
   incluye la relajación que faltaba y no significa nada.

4. **El orden de los átomos.** Las dos estructuras tienen que listar los
   átomos en el MISMO orden: la imagen i se interpola átomo por átomo. Si
   el orden no coincide, la interpolación hace que los átomos se
   atraviesen. Olla-DFT lo comprueba y se niega.

LO QUE NO ES UNA BARRERA
------------------------
La diferencia de energía que sale de aquí es electrónica, a 0 K y sin
punto cero. La barrera que compara con una constante de velocidad
experimental necesita la corrección de punto cero y la entropía — eso
está en `olla-dft thermochem`.
"""

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from qekit.core import provenance
from qekit.core import style as qstyle
from qekit.core.errors import ErrorDeUso

BOHR_A = 0.529177210903


@dataclass
class NebRun:
    s: np.ndarray = None            # coordenada de reacción (normalizada)
    energias: np.ndarray = None     # eV, relativas a la primera imagen
    fuerzas: np.ndarray = None      # eV/Å, proyectadas sobre el camino
    s_int: np.ndarray = None        # interpolación cúbica
    e_int: np.ndarray = None
    barrera_ida: float = None       # eV
    barrera_vuelta: float = None
    delta_E: float = None           # producto - reactivo
    imagen_cima: int = None
    n_imagenes: int = 0
    ci: bool = None
    convergido: bool = None
    pasos: int = None
    imagenes_malas: list = field(default_factory=list)
    avisos: list = field(default_factory=list)


# ----------------------------------------------------------------------
# Preparación
# ----------------------------------------------------------------------
def comprobar_extremos(inicial, final) -> list:
    """El error que arruina un NEB antes de empezar: extremos incompatibles."""
    problemas = []
    s1 = list(inicial.get_chemical_symbols())
    s2 = list(final.get_chemical_symbols())
    if len(s1) != len(s2):
        problemas.append(
            f"the two structures have a different number of atoms "
            f"({len(s1)} and {len(s2)}). A reaction path conserves the "
            "atoms.")
        return problemas
    if s1 != s2:
        distintos = [i for i, (a, b) in enumerate(zip(s1, s2)) if a != b]
        problemas.append(
            f"the atoms are not in the same ORDER: they differ at "
            f"positions {distintos[:6]}"
            f"{'...' if len(distintos) > 6 else ''}.\n"
            "The interpolation goes atom by atom, so with the order "
            "changed the atoms pass\nthrough each other and the path means "
            "nothing. Reorder one of the two.")
    c1 = np.array(inicial.get_cell())
    c2 = np.array(final.get_cell())
    if not np.allclose(c1, c2, atol=1e-4):
        problemas.append(
            "the cells are not equal. neb.x does not relax the cell: both "
            "endpoints must\nbe in the same one.")
    d = np.linalg.norm(final.get_positions() - inicial.get_positions(), axis=1)
    if d.max() < 1e-6:
        problemas.append(
            "the two structures are identical: there is no path to search for.")
    return problemas


def build_neb_input(inicial, final, pseudos: dict, prefix: str,
                    pseudo_dir: str, ecutwfc: float, ecutrho: float,
                    kcard: str, n_imagenes: int = 7, ci: bool = True,
                    path_thr: float = 0.05, nstep_path: int = 50,
                    opt_scheme: str = "broyden", k_min: float = 0.2,
                    k_max: float = 0.3, insulator: bool = True,
                    degauss: float = 0.01, smearing: str = "cold",
                    nspin: int = 1, magnetization: dict = None,
                    intermedias=None, fijos=None) -> str:
    """Arma el input de neb.x, con su estructura de BEGIN/END."""
    from qekit.modules import inputgen

    motor = inputgen.build_pw_input(
        atoms=inicial, pseudos=pseudos, calculation="scf", prefix=prefix,
        pseudo_dir=pseudo_dir, ecutwfc=ecutwfc, ecutrho=ecutrho,
        kcard=kcard, insulator=insulator, degauss=degauss,
        smearing=smearing, nspin=nspin, magnetization=magnetization)
    # neb.x lleva su propio bloque de posiciones por imagen: se quitan del
    # input del motor las tarjetas de estructura y se dejan solo las
    # namelists y ATOMIC_SPECIES.
    motor = _recortar_motor(motor)

    lineas = ["BEGIN", "BEGIN_PATH_INPUT", "&PATH",
              "  restart_mode      = 'from_scratch',",
              "  string_method     = 'neb',",
              f"  nstep_path        = {nstep_path},",
              "  ds                = 1.D0,",
              f"  opt_scheme        = '{opt_scheme}',",
              f"  num_of_images     = {n_imagenes},",
              f"  k_max             = {k_max}D0,",
              f"  k_min             = {k_min}D0,",
              f"  CI_scheme         = '{'auto' if ci else 'no-CI'}',",
              f"  path_thr          = {path_thr}D0,",
              "/", "END_PATH_INPUT",
              "BEGIN_ENGINE_INPUT", motor.rstrip(),
              "BEGIN_POSITIONS", "FIRST_IMAGE",
              _tarjeta_posiciones(inicial, fijos)]
    for extra in (intermedias or []):
        lineas += ["INTERMEDIATE_IMAGE", _tarjeta_posiciones(extra, fijos)]
    lineas += ["LAST_IMAGE", _tarjeta_posiciones(final, fijos),
               "END_POSITIONS",
               _tarjeta_celda(inicial),
               "END_ENGINE_INPUT", "END"]
    return "\n".join(lineas) + "\n"


def _recortar_motor(texto: str) -> str:
    """Deja las namelists, ATOMIC_SPECIES y K_POINTS; quita las posiciones."""
    fuera, saltando = [], False
    for ln in texto.split("\n"):
        if ln.startswith("ATOMIC_POSITIONS") or ln.startswith("CELL_PARAMETERS"):
            saltando = True
            continue
        if ln.startswith(("K_POINTS", "ATOMIC_SPECIES")):
            saltando = False
        if not saltando:
            fuera.append(ln)
    return "\n".join(fuera)


def _tarjeta_posiciones(atoms, fijos=None) -> str:
    lineas = ["ATOMIC_POSITIONS { angstrom }"]
    fijos = set(fijos or [])
    for i, (s, p) in enumerate(zip(atoms.get_chemical_symbols(),
                                   atoms.get_positions())):
        cong = "  0 0 0" if i in fijos else ""
        lineas.append(f"  {s:3s} {p[0]:14.9f} {p[1]:14.9f} {p[2]:14.9f}{cong}")
    return "\n".join(lineas)


def _tarjeta_celda(atoms) -> str:
    c = np.array(atoms.get_cell())
    lineas = ["CELL_PARAMETERS { angstrom }"]
    for v in c:
        lineas.append(f"  {v[0]:14.9f} {v[1]:14.9f} {v[2]:14.9f}")
    return "\n".join(lineas)


def prepare(inicial, final, outdir: str = "neb", n_imagenes: int = 7,
            ci: bool = True, pseudo_dir: str = None, ecutwfc: float = None,
            ecutrho: float = None, kspacing: float = None,
            insulator: bool = True, path_thr: float = 0.05,
            nstep_path: int = 50, nspin: int = 1, magnetization: dict = None,
            fijos=None) -> tuple:
    from qekit.modules import sweep

    problemas = comprobar_extremos(inicial, final)
    if problemas:
        raise ErrorDeUso("the endpoints of the path are not compatible:\n\n" +
                         "\n\n".join("  " + p for p in problemas))

    common = sweep.prepare_common(inicial, pseudo_dir, ecutwfc, ecutrho,
                                  insulator)
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    grid = sweep.default_grid(inicial, kspacing)
    texto = build_neb_input(
        inicial, final, common["pseudos"], common["prefix"],
        common["pseudo_dir"], common["ecutwfc"], common["ecutrho"],
        kcard=f"K_POINTS automatic\n  {grid[0]} {grid[1]} {grid[2]} 0 0 0\n",
        n_imagenes=n_imagenes, ci=ci, path_thr=path_thr,
        nstep_path=nstep_path, insulator=insulator,
        degauss=common["degauss"], smearing=common["smearing"],
        nspin=nspin, magnetization=magnetization, fijos=fijos)
    sweep.write_input(out / "neb.in", texto)

    d = np.linalg.norm(final.get_positions() - inicial.get_positions(), axis=1)
    mueve = int(np.sum(d > 0.1))
    rep = ["--- Reaction path (NEB) ---",
           f"Structure: {inicial.get_chemical_formula()} "
           f"({len(inicial)} atoms)",
           f"Atoms moving more than 0.1 Å: {mueve}  "
           f"(largest displacement: {d.max():.2f} Å)",
           f"Images: {n_imagenes}   climbing image: "
           f"{'yes' if ci else 'NO'}",
           f"Path force threshold: {path_thr} eV/Å",
           "",
           f"File in '{out.resolve()}': neb.in",
           "Run it with:  neb.x -inp neb.in > neb.out",
           ""]
    if not ci:
        rep += ["WARNING: without a climbing image the barrier is UNDERESTIMATED. The "
                "highest image\nstays close to the maximum but not on top of it, "
                "and the typical error is\ntenths of an eV.", ""]
    if n_imagenes < 5:
        rep += [f"WARNING: {n_imagenes} images is too few. With so few the "
                "maximum can fall between\ntwo images and be missed.", ""]
    rep += ["Both endpoints must be RELAXED with these same "
            "parameters.\nOtherwise, the barrier includes the missing relaxation."]
    warn = sweep.missing_pseudo_warning(common)
    if warn:
        rep.append(warn)
    return common, "\n".join(rep)


# ----------------------------------------------------------------------
# Lectura
# ----------------------------------------------------------------------
def collect(path, prefix: str = None) -> NebRun:
    p = Path(path)
    dats = sorted(p.glob("*.dat"))
    dats = [d for d in dats if not d.name.endswith("_int.dat")]
    if prefix:
        dats = [d for d in dats if d.stem == prefix] or dats
    if not dats:
        raise ErrorDeUso(
            f"there is no <prefix>.dat in {p}. neb.x writes it with the "
            "energy profile;\nif it is not there, check neb.out.")
    d = np.loadtxt(dats[0])
    if d.ndim == 1:
        d = d.reshape(1, -1)
    run = NebRun(s=d[:, 0], energias=d[:, 1],
                 fuerzas=d[:, 2] if d.shape[1] > 2 else None,
                 n_imagenes=len(d))

    intp = dats[0].with_suffix(".int")
    if intp.exists():
        di = np.loadtxt(intp)
        if di.ndim == 2 and di.shape[1] >= 2:
            run.s_int, run.e_int = di[:, 0], di[:, 1]

    salida = sorted(p.glob("*.out"))
    if salida:
        texto = salida[0].read_text(errors="ignore")
        m = re.search(r"activation energy \(->\)\s*=\s*(-?[\d.]+)\s*eV", texto)
        if m:
            run.barrera_ida = float(m.group(1))
        m = re.search(r"activation energy \(<-\)\s*=\s*(-?[\d.]+)\s*eV", texto)
        if m:
            run.barrera_vuelta = float(m.group(1))
        run.convergido = "neb: convergence achieved" in texto.lower() or \
            "path length" in texto and "convergence achieved" in texto.lower()
        pasos = re.findall(r"iteration\s*=\s*(\d+)", texto)
        if pasos:
            run.pasos = int(pasos[-1])
        m = re.search(r"CI_scheme\s*=\s*(\S+)", texto)
        if m:
            run.ci = m.group(1).strip().lower() != "no-ci"
        # Un scf que no converge en una imagen envenena su energia y la
        # fuerza, y neb.x sigue adelante con un WARNING que se pierde entre
        # cientos de lineas. Es la causa mas comun de un perfil dentado.
        malas = sorted(set(int(x) for x in re.findall(
            r"scf convergence NOT achieved on image\s+(\d+)", texto)))
        if malas:
            run.imagenes_malas = malas
            run.avisos.append(
                "The scf did NOT converge on image(s) "
                + ", ".join(str(i) for i in malas) +
                ".\nThe energy and force of those images are not reliable, and "
                "the whole path\ninherits the error: that is why the profile comes out "
                "jagged. Lower mixing_beta, raise\nelectron_maxstep, or loosen "
                "the engine's conv_thr.")

    e = run.energias
    run.imagen_cima = int(np.argmax(e))
    if run.barrera_ida is None:
        run.barrera_ida = float(e.max() - e[0])
    if run.barrera_vuelta is None:
        run.barrera_vuelta = float(e.max() - e[-1])
    run.delta_E = float(e[-1] - e[0])
    return run


def report(run: NebRun) -> str:
    lines = ["--- Reaction path (NEB) ---",
             f"Images: {run.n_imagenes}"]
    if run.ci is not None:
        lines.append(f"Climbing image: {'yes' if run.ci else 'NO'}")
    if run.pasos is not None:
        lines.append(f"Path iterations: {run.pasos}")
    if run.convergido is not None:
        lines.append(f"Converged: {'yes' if run.convergido else 'NO'}")
    lines += ["",
              f"Forward barrier   (reactant -> product): "
              f"{run.barrera_ida:8.4f} eV  "
              f"({run.barrera_ida * 96.485:7.1f} kJ/mol)",
              f"Reverse barrier   (product -> reactant): "
              f"{run.barrera_vuelta:8.4f} eV  "
              f"({run.barrera_vuelta * 96.485:7.1f} kJ/mol)",
              f"Reaction energy (product - reactant): "
              f"{run.delta_E:+8.4f} eV  "
              f"({run.delta_E * 96.485:+7.1f} kJ/mol)",
              "",
              f"{'image':>7s} {'s':>8s} {'E (eV)':>10s} "
              f"{'F (eV/Å)':>10s}"]
    for i in range(run.n_imagenes):
        f = run.fuerzas[i] if run.fuerzas is not None else float("nan")
        marca = "  <- top" if i == run.imagen_cima else ""
        lines.append(f"{i + 1:7d} {run.s[i]:8.4f} {run.energias[i]:10.4f} "
                     f"{f:10.4f}{marca}")

    # ¿el máximo interpolado cae lejos de toda imagen calculada?
    if run.s_int is not None and run.e_int is not None:
        smax = float(run.s_int[int(np.argmax(run.e_int))])
        cerca = float(np.min(np.abs(run.s - smax)))
        paso = float(np.mean(np.diff(run.s))) if run.n_imagenes > 1 else 1.0
        if cerca > 0.4 * paso:
            lines += ["",
                      "The maximum of the interpolated curve falls between two "
                      "computed images.\nThat means the transition "
                      "state is NOT sampled: raise\n"
                      "--images, or let the climbing image settle "
                      "on top of it."]
    if run.ci is False:
        lines += ["",
                  "Without a climbing image, this barrier is a LOWER BOUND: "
                  "the highest image\nstays below the real maximum."]
    for a in run.avisos:
        lines += ["", a]
    if run.convergido is False:
        lines += ["",
                  "The path did NOT converge. The numbers above are "
                  "provisional: raise\nnstep_path, or loosen path_thr and "
                  "tighten it again in a second stage."]
    lines += ["",
              "This barrier is ELECTRONIC, at 0 K and without zero-point "
              "energy. To\ncompare it with an experimental activation "
              "energy, the thermal\ncorrections are needed: "
              "'olla-dft thermochem'."]
    return "\n".join(lines)


def export(run: NebRun, outdir: str = ".") -> list:
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    f = out / "NEB.dat"
    cols = [run.s, run.energias]
    nombres = ["s", "E(eV)"]
    if run.fuerzas is not None:
        cols.append(run.fuerzas); nombres.append("F(eV/A)")
    np.savetxt(f, np.column_stack(cols), fmt="%14.6f",
               header=provenance.header_plain(
                   "reaction path",
                   {"barrera_ida_eV": round(run.barrera_ida, 5),
                    "barrera_vuelta_eV": round(run.barrera_vuelta, 5),
                    "delta_E_eV": round(run.delta_E, 5),
                    "imagenes": run.n_imagenes, "CI": run.ci},
                   titulo="NEB energy profile") + "\n" +
               "  ".join(f"{n:>14s}" for n in nombres), comments="# ")
    txt = out / "NEB.txt"
    txt.write_text(report(run) + "\n", encoding="utf-8")
    return [str(f), str(txt)]


def plot(run: NebRun, outfile: str = "neb", formats="pdf,png",
         theme: str = None, family: str = None, background: str = None,
         palette=None, usetex: bool = None, width="single",
         journal: str = "generic", mono: bool = False,
         dpi: int = None) -> list:
    try:
        import matplotlib
        matplotlib.use("Agg")
    except ImportError as exc:                          # pragma: no cover
        raise RuntimeError("matplotlib is not installed.") from exc

    st = qstyle.apply(theme, family=family, background=background,
                      palette=palette, usetex=usetex, mono=mono)
    fig, ax = qstyle.new_figure(width, journal, 0.72)
    colores = qstyle.palette(4, mono=mono)

    if run.s_int is not None:
        ax.plot(run.s_int, run.e_int, lw=1.4, color=colores[0],
                label="interpolation")
    ax.plot(run.s, run.energias, "o", ms=4.5, color=colores[1],
            label="images", zorder=3)
    ax.plot(run.s[run.imagen_cima], run.energias[run.imagen_cima], "o",
            ms=7, mfc="none", mew=1.4, color=colores[2],
            label="transition state", zorder=4)
    ax.axhline(0.0, color=qstyle.INK_FAINT, lw=st["axis_line"],
               dashes=[3.5, 2.0])
    ax.annotate("", xy=(run.s[run.imagen_cima], run.energias.max()),
                xytext=(run.s[run.imagen_cima], run.energias[0]),
                arrowprops=dict(arrowstyle="<->", lw=0.8,
                                color=qstyle.INK_FAINT))
    ax.text(run.s[run.imagen_cima], 0.5 * (run.energias.max() +
                                           run.energias[0]),
            f"  {run.barrera_ida:.3f} eV", va="center", fontsize="small")
    ax.set_xlabel("reaction coordinate")
    ax.set_ylabel("E (eV)")
    ax.legend(frameon=False, fontsize="small")
    return qstyle.save(fig, outfile, formats, dpi=dpi)
