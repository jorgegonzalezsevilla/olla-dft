# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Sitios de adsorción sobre una losa: enumerarlos, montarlos y compararlos.

La energía de adsorción es una resta de tres números:

    E_ads = E(losa + molécula) - E(losa) - n·E(molécula)

La resta ya vivía en `thermochem.adsorcion()`. Lo que faltaba, y es lo que
de verdad cuesta a mano, es todo lo de alrededor: encontrar los sitios de
la superficie, no repetir los que la simetría hace equivalentes, colocar la
molécula a una altura sensata en cada uno, y garantizar que los tres
cálculos son comparables. Ese último punto es el que más silenciosamente se
rompe: si la losa limpia y la losa con molécula no comparten celda, cutoff,
malla k y pseudos, la resta da un número perfectamente formado que no
significa nada. Aquí los tres cálculos se generan a la vez y de la misma
plantilla, precisamente para que no puedan divergir.

Tipos de sitio que se enumeran:
  top      encima de un átomo de la superficie
  bridge   sobre el punto medio de dos átomos vecinos
  hollow   sobre el centro de un triángulo de átomos (fcc, hcp, 4-fold...)
"""

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from qekit.core import provenance, qeout, structure
from qekit.core.errors import ErrorDeUso, FaltanDatos
from qekit.core import style as qstyle
from qekit.modules import sweep, thermochem

TIPOS = ("top", "bridge", "hollow")

# Distancia máxima (Å) para considerar vecinos dos átomos de la superficie al
# formar los puentes. Por encima de esto el "puente" no está entre nada.
R_VECINO = 3.6

# Tolerancia (Å) para decidir qué átomos forman la capa superior.
TOL_CAPA = 0.6


@dataclass
class Sitio:
    tipo: str
    xy: tuple                   # posición cartesiana en el plano (Å)
    z: float                    # altura de la superficie bajo el sitio (Å)
    huella: object = None       # firma de distancias, para deduplicar
    rotacion: float = 0.0       # grados alrededor de la normal
    etiqueta: str = ""


@dataclass
class AdsorbRun:
    jobs: list = field(default_factory=list)
    sitios: list = field(default_factory=list)
    energies: list = field(default_factory=list)     # eV, uno por sitio
    converged: list = field(default_factory=list)
    alturas: list = field(default_factory=list)      # Å tras relajar
    contactos: list = field(default_factory=list)    # Å, distancia mínima
    E_slab: float = None
    E_mol: float = None
    slab_ok: bool = None
    mol_ok: bool = None
    molecula: str = ""
    n_mol: int = 1
    natoms_slab: int = 0
    natoms_mol: int = 1
    altura_inicial: float = None
    relajado: bool = True

    @property
    def energias_ads(self) -> list:
        if self.E_slab is None or self.E_mol is None:
            return [None] * len(self.sitios)
        out = []
        for e in self.energies:
            out.append(None if e is None
                       else thermochem.adsorcion(e, self.E_slab, self.E_mol,
                                                 n=self.n_mol)["E_ads"])
        return out


# ----------------------------------------------------------------------
# Geometría de la superficie
# ----------------------------------------------------------------------
def atomos_superficie(slab, cara: str = "top", tol: float = TOL_CAPA) -> list:
    """Índices de los átomos de la capa expuesta."""
    z = slab.get_positions()[:, 2]
    if cara == "top":
        ref = z.max()
        return [i for i in range(len(slab)) if z[i] >= ref - tol]
    if cara == "bottom":
        ref = z.min()
        return [i for i in range(len(slab)) if z[i] <= ref + tol]
    raise ErrorDeUso(f"--face must be 'top' or 'bottom'; got '{cara}'.")


def _replicas(slab, idx, n=1):
    """Posiciones de los átomos `idx` y sus réplicas periódicas en el plano."""
    cell = slab.cell.array
    pos = slab.get_positions()[idx]
    out, orig = [], []
    for i in range(-n, n + 1):
        for j in range(-n, n + 1):
            desp = i * cell[0] + j * cell[1]
            for k, p in enumerate(pos):
                out.append(p + desp)
                orig.append(idx[k])
    return np.array(out), np.array(orig)


# Cuántas distancias del entorno se comparan para decidir si dos sitios son
# el mismo, y con qué tolerancia en Å. Un número FIJO de vecinos y no un
# radio: con un radio, un átomo que cae justo en el borde entra en la firma
# de un sitio y no en la de su equivalente, y dos sitios idénticos salen
# distintos. Se comprobó en Al(111) 2x2, donde el radio daba 3 huecos
# distintos donde solo hay dos (fcc y hcp).
N_VECINOS_HUELLA = 24
TOL_HUELLA = 0.05


def _huella(slab, punto, k: int = N_VECINOS_HUELLA) -> np.ndarray:
    """Firma de un sitio: distancias a sus k vecinos más cercanos, ordenadas.

    Dos sitios equivalentes por simetría tienen el mismo entorno. Se incluyen
    TODAS las capas, no solo la superficial, porque es lo único que separa un
    hueco fcc de uno hcp: por arriba se ven idénticos y difieren en si hay
    un átomo debajo en la segunda capa.
    """
    cell = slab.cell.array
    anchos = [np.linalg.norm(cell[0]), np.linalg.norm(cell[1])]
    # réplicas suficientes para que la esfera de vecinos quepa entera
    n = max(2, int(np.ceil(12.0 / max(1e-6, min(anchos)))) + 1)
    todos, _ = _replicas(slab, list(range(len(slab))), n=n)
    d = np.sort(np.linalg.norm(todos - np.asarray(punto), axis=1))
    if len(d) < k:
        d = np.pad(d, (0, k - len(d)), constant_values=d[-1] if len(d) else 0.0)
    return d[:k]


def _misma_huella(a: np.ndarray, b: np.ndarray, tol: float = TOL_HUELLA) -> bool:
    return bool(np.max(np.abs(a - b)) < tol)


def sitios(slab, cara: str = "top", tipos=TIPOS, tol: float = TOL_CAPA,
           r_vecino: float = R_VECINO) -> list:
    """Enumera los sitios de adsorción no equivalentes de la cara expuesta."""
    from scipy.spatial import Delaunay

    idx = atomos_superficie(slab, cara, tol)
    if not idx:
        raise ErrorDeUso("no surface atoms found; is this a slab with "
                         "vacuum? Cut one with 'olla-dft surface'.")
    pos = slab.get_positions()[idx]
    z_sup = pos[:, 2].max() if cara == "top" else pos[:, 2].min()
    rep, orig = _replicas(slab, idx, n=1)

    cand = []
    if "top" in tipos:
        for p in pos:
            cand.append(Sitio("top", (p[0], p[1]), p[2]))

    if "bridge" in tipos:
        for a in range(len(pos)):
            d = np.linalg.norm(rep - pos[a], axis=1)
            for b in np.where((d > 0.1) & (d <= r_vecino))[0]:
                m = (pos[a] + rep[b]) / 2.0
                cand.append(Sitio("bridge", (m[0], m[1]), float(m[2])))

    if "hollow" in tipos:
        plano = rep[:, :2]
        if len(plano) >= 4:
            try:
                tri = Delaunay(plano)
            except Exception:                               # noqa: BLE001
                tri = None
            if tri is not None:
                for simplex in tri.simplices:
                    v = rep[simplex]
                    lados = [np.linalg.norm(v[i] - v[j])
                             for i, j in ((0, 1), (1, 2), (0, 2))]
                    if max(lados) > r_vecino * 1.6:
                        continue    # triángulo estirado: no es un hueco real
                    c = v.mean(axis=0)
                    cand.append(Sitio("hollow", (c[0], c[1]), float(c[2])))

    # --- quedarse solo con los que están dentro de la celda y no se repiten ---
    cell2 = slab.cell.array[:2, :2]
    inv = np.linalg.inv(cell2.T)
    unicos = []
    for s in cand:
        f = inv @ np.array(s.xy)
        if not (-1e-6 <= f[0] < 1 - 1e-6 and -1e-6 <= f[1] < 1 - 1e-6):
            # traer a la celda de referencia antes de comparar
            f = f % 1.0
            xy = cell2.T @ f
            s = Sitio(s.tipo, (float(xy[0]), float(xy[1])), s.z)
        h = _huella(slab, (s.xy[0], s.xy[1], z_sup))
        if any(_misma_huella(h, u.huella) for u in unicos):
            continue
        s.huella = h
        unicos.append(s)

    orden = {"top": 0, "bridge": 1, "hollow": 2}
    unicos.sort(key=lambda s: (orden[s.tipo], s.xy))
    cuenta = {}
    for s in unicos:
        cuenta[s.tipo] = cuenta.get(s.tipo, 0) + 1
        s.etiqueta = f"{s.tipo}{cuenta[s.tipo]}"
    return unicos


def cargar_molecula(nombre: str):
    """Molécula por nombre (base de ASE) o desde un archivo."""
    p = Path(nombre)
    if p.exists():
        return structure.load(str(p))
    try:
        from ase.build import molecule
        return molecule(nombre)
    except Exception:                                       # noqa: BLE001
        try:
            from ase.collections import g2
            disponibles = ", ".join(sorted(g2.names)[:14])
        except Exception:                                   # noqa: BLE001
            disponibles = "CO, CO2, H2O, NH3, O2, CH4..."
        raise ErrorDeUso(
            f"'{nombre}' is neither an existing file nor a molecule in the "
            f"ASE database. Some that are: {disponibles}. You can also "
            f"pass an .xyz or .cif file with the molecule.") from None


def colocar(slab, mol, sitio: Sitio, altura: float = 2.0,
            rotacion: float = 0.0, ancla: int = 0, cara: str = "top"):
    """Devuelve la losa con la molécula puesta en el sitio."""
    m = mol.copy()
    m.set_cell(slab.cell)
    m.set_pbc(slab.pbc)
    if rotacion:
        m.rotate(rotacion, "z", center=m.get_positions()[ancla])
    p = m.get_positions()
    ancla_pos = p[ancla]
    signo = 1.0 if cara == "top" else -1.0
    destino = np.array([sitio.xy[0], sitio.xy[1], sitio.z + signo * altura])
    m.set_positions(p + (destino - ancla_pos))
    out = slab.copy()
    out += m
    return out


# ----------------------------------------------------------------------
# Preparación
# ----------------------------------------------------------------------
def prepare(slab, molecula: str, outdir: str = "adsorb",
            altura: float = 2.0, tipos=TIPOS, cara: str = "top",
            rotaciones: int = 1, ancla: int = 0, pseudo_dir: str = None,
            insulator: bool = False, ecutwfc: float = None,
            ecutrho: float = None, kspacing: float = None,
            relax_ions: bool = True, vdw: str = None, dipolo: bool = False,
            nspin: int = 1, magnetization: dict = None) -> tuple:
    from qekit.core import kpoints as kp

    if 2 not in kp.direcciones_con_vacio(slab):
        raise ErrorDeUso(
            "this structure has no vacuum along c: the adsorption energy "
            "needs a slab with vacuum above it. Cut one with "
            "'olla-dft surface -m \"1 1 1\" --vacuum 20'.")
    if rotaciones < 1:
        raise ErrorDeUso(f"--rotations must be at least 1; got {rotaciones}.")
    tipos = tuple(tipos)
    malos = [t for t in tipos if t not in TIPOS]
    if malos:
        raise ErrorDeUso(
            f"unknown site type: {', '.join(malos)}. "
            f"Options: {', '.join(TIPOS)}.")

    mol = cargar_molecula(molecula)
    if ancla >= len(mol):
        raise ErrorDeUso(
            f"--anchor {ancla} does not exist: the molecule has {len(mol)} atoms "
            f"(numbered from 0).")

    lista = sitios(slab, cara=cara, tipos=tipos)
    if rotaciones > 1 and len(mol) > 1:
        ampliada = []
        for s in lista:
            for k in range(rotaciones):
                s2 = Sitio(s.tipo, s.xy, s.z, s.huella,
                           rotacion=360.0 * k / rotaciones,
                           etiqueta=f"{s.etiqueta}_r{k}")
                ampliada.append(s2)
        lista = ampliada
    elif rotaciones > 1:
        # una molécula de un átomo no cambia al girar; girarla sería pagar
        # N veces el mismo cálculo para obtener N veces el mismo número.
        rotaciones = 1

    # Los pseudos y los cutoffs se resuelven sobre la UNIÓN de losa y
    # molécula: si se hicieran solo con la losa, el carbono y el oxígeno del
    # adsorbato no tendrían pseudo, y el cutoff sería el de los metales, que
    # es demasiado bajo para el oxígeno. Los tres cálculos comparten el
    # resultado, que es justo lo que hace la resta comparable.
    conjunto = slab + mol
    conjunto.set_cell(slab.cell)
    common = sweep.prepare_common(conjunto, pseudo_dir, ecutwfc, ecutrho,
                                  insulator,
                                  prefix=slab.get_chemical_formula(
                                      mode="hill", empirical=True))
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    grid = sweep.default_grid(slab, kspacing)
    calc = "relax" if relax_ions else "scf"

    run = AdsorbRun(sitios=lista, molecula=molecula, natoms_slab=len(slab),
                    natoms_mol=len(mol), altura_inicial=altura,
                    relajado=relax_ions)

    # Los tres tipos de cálculo comparten celda, cutoffs, malla k y pseudos:
    # es la única forma de que la resta signifique algo. La molécula va en la
    # MISMA celda que la losa, no en una caja aparte, por lo mismo.
    mol_sola = mol.copy()
    mol_sola.set_cell(slab.cell)
    mol_sola.set_pbc(slab.pbc)
    mol_sola.center()

    # --dipole: la sierra de la corrección dipolar (tefield/dipfield, edir=3)
    # va en los TRES cálculos, no solo en el de la losa con adsorbato: si la
    # referencia se calcula sin corregir, la resta arrastra el error.
    extras = dict(vdw=vdw, nspin=nspin, magnetization=magnetization,
                  dipole_correction=3 if dipolo else False)
    run.jobs.append(sweep.write_scf_job(
        slab, common, out / "_losa", "clean slab", grid,
        meta={"papel": "slab"}, calculation=calc, **extras))
    run.jobs.append(sweep.write_scf_job(
        mol_sola, common, out / "_molecula", f"isolated {molecula}", grid,
        meta={"papel": "mol"}, calculation=calc, **extras))

    for s in lista:
        sistema = colocar(slab, mol, s, altura=altura, rotacion=s.rotacion,
                          ancla=ancla, cara=cara)
        job = sweep.write_scf_job(
            sistema, common, out / s.etiqueta, s.etiqueta, grid,
            meta={"papel": "ads", "sitio": s.etiqueta}, calculation=calc,
            **extras)
        run.jobs.append(job)

    sweep.write_run_script(run.jobs, out / "run.sh")

    cuenta = {}
    for s in lista:
        cuenta[s.tipo] = cuenta.get(s.tipo, 0) + 1
    report = ["--- Adsorption sites ---",
              f"Slab: {slab.get_chemical_formula()} ({len(slab)} atoms), "
              f"face {cara}",
              f"Adsorbate: {molecula} ({len(mol)} atoms), "
              f"anchor = atom {ancla} ({mol.get_chemical_symbols()[ancla]}), "
              f"initial height {altura:g} Å",
              "Non-equivalent sites: "
              + ", ".join(f"{n} {t}" for t, n in sorted(cuenta.items()))
              + f"  ({len(lista)} calculations"
              + (f", {rotaciones} rotations each" if rotaciones > 1 else "")
              + ")",
              f"k-mesh: {grid[0]}x{grid[1]}x{grid[2]}  |  "
              + ("relaxed positions" if relax_ions else "fixed positions"),
              "References: clean slab and isolated molecule, in the SAME cell "
              "and with the\n  same cutoffs and mesh, so that the difference is "
              "valid."]
    if not vdw:
        report.append(
            "WARNING: no van der Waals correction. In physisorption (closed-shell "
            "molecules\n  on surfaces) the bond IS dispersion: without --vdw "
            "the energy comes out\n  near zero and the geometry unbound.")
    if not dipolo and cara == "top":
        report.append(
            "Suggestion: a molecule adsorbed on a single face leaves the slab "
            "polar.\n  With --dipole the artificial dipole across the vacuum "
            "is cancelled.")
    warn = sweep.missing_pseudo_warning(common)
    if warn:
        report.append(warn)
    report += ["", f"{len(run.jobs)} calculations written to '{out.resolve()}'",
               "Run them with --run, or by hand with ./run.sh inside that folder."]
    return run, "\n".join(report)


# ----------------------------------------------------------------------
# Recolección
# ----------------------------------------------------------------------
def _leer(job, por_dir):
    r = por_dir.get(str(job.directory))
    if r is not None and r.ok and r.result is not None:
        return r.result
    try:
        return qeout.read_xml(str(job.directory))
    except Exception:                                       # noqa: BLE001
        return None


def collect(run: AdsorbRun, results: list = None) -> AdsorbRun:
    por_dir = {str(r.job.directory): r for r in (results or [])}
    run.energies, run.converged, run.alturas, run.contactos = [], [], [], []
    n_slab = run.natoms_slab
    for job in run.jobs:
        papel = job.meta.get("papel")
        res = _leer(job, por_dir)
        if papel == "slab":
            run.E_slab = res.total_energy if res else None
            run.slab_ok = res.converged if res else None
            continue
        if papel == "mol":
            run.E_mol = res.total_energy if res else None
            run.mol_ok = res.converged if res else None
            continue
        run.energies.append(res.total_energy if res else None)
        run.converged.append(res.converged if res else None)
        if res is not None and res.positions is not None and len(res.positions) > n_slab:
            p = np.asarray(res.positions)
            losa, ads = p[:n_slab], p[n_slab:]
            run.alturas.append(float(ads[:, 2].min() - losa[:, 2].max()))
            d = np.linalg.norm(ads[:, None, :] - losa[None, :, :], axis=2)
            run.contactos.append(float(d.min()))
        else:
            run.alturas.append(None)
            run.contactos.append(None)
    return run


# ----------------------------------------------------------------------
# Reporte
# ----------------------------------------------------------------------
def report(run: AdsorbRun) -> str:
    if not run.energies:
        raise FaltanDatos(
            "no results yet. Run the calculations (--run, or ./run.sh "
            "in the folder) and come back with --collect.")
    L = ["--- Adsorption energies ---",
         f"Adsorbate: {run.molecula}   |   slab of {run.natoms_slab} atoms"]
    faltan = []
    if run.E_slab is None:
        faltan.append("the clean slab")
    if run.E_mol is None:
        faltan.append("the isolated molecule")
    if faltan:
        L.append("")
        L.append("E_ads cannot be computed: missing the energy of "
                 + " and ".join(faltan) + ".")
        L.append("  Without both references there is no difference to take; run "
                 "those two calculations\n  (they are in _losa and _molecula) and "
                 "come back with --collect.")
        return "\n".join(L)

    L.append(f"E(slab) = {run.E_slab:.6f} eV    "
             f"E({run.molecula}) = {run.E_mol:.6f} eV")
    L.append("")
    L.append(f"  {'site':<12s} {'E_ads (eV)':>11s} {'height (Å)':>11s} "
             f"{'contact (Å)':>13s}")
    L.append("  " + "-" * 51)

    eads = run.energias_ads
    filas = sorted(
        [(i, s) for i, s in enumerate(run.sitios)],
        key=lambda t: (eads[t[0]] is None, eads[t[0]] if eads[t[0]] is not None else 0.0))
    for i, s in filas:
        e = eads[i]
        if e is None:
            L.append(f"  {s.etiqueta:<12s} {'no result':>11s}")
            continue
        alt = run.alturas[i]
        con = run.contactos[i]
        fila = (f"  {s.etiqueta:<12s} {e:>11.4f} "
                f"{(alt if alt is not None else float('nan')):>11.3f} "
                f"{(con if con is not None else float('nan')):>13.3f}")
        if run.converged[i] is False:
            fila += "   << NOT CONVERGED"
        L.append(fila)

    validos = [(i, eads[i]) for i, _ in filas if eads[i] is not None]
    if validos:
        mejor_i, mejor_e = validos[0]
        s = run.sitios[mejor_i]
        L.append("")
        L.append(f"Most favourable site: {s.etiqueta} ({s.tipo}), "
                 f"E_ads = {mejor_e:.4f} eV")
        if mejor_e > 0:
            cerca = [run.contactos[i] for i, _ in filas
                     if run.contactos[i] is not None]
            if not run.relajado:
                L.append(
                    "  POSITIVE at every site, and the calculation used "
                    "FIXED positions:\n  most likely the initial height "
                    f"({run.altura_inicial:g} Å) is not the equilibrium "
                    "one\n  and you are measuring the repulsion. Remove "
                    "--fixed-ions to let it relax.")
            elif cerca and min(cerca) < 1.2:
                L.append(
                    f"  POSITIVE, and the shortest contact is {min(cerca):.2f} Å: "
                    "the atoms are\n  sitting on top of each other. Check the "
                    "geometry before concluding anything.")
            else:
                L.append("  POSITIVE: at this level of theory the adsorbate does "
                         "NOT bind at any\n  site tried. If you expected "
                         "physisorption, try --vdw.")
        elif mejor_e > -0.30:
            L.append("  Weak physisorption (|E_ads| < 0.3 eV): at room "
                     "temperature the\n  molecule desorbs. The number depends "
                     "strongly on the dispersion correction.")
        elif mejor_e < -2.0 and run.natoms_mol > 1:
            L.append("  Very strong bond (|E_ads| > 2 eV): usually there is a "
                     "reaction, not\n  molecular adsorption. Look at the relaxed "
                     "geometry: the molecule may have\n  dissociated and you are "
                     "measuring the energy of the fragments.")
        elif mejor_e < -2.0:
            L.append("  Strong chemisorption. Mind the reference: this "
                     "E_ads is measured against\n  the isolated ATOM, not the "
                     "molecule. To compare with the\n  literature on "
                     "diatomic molecules you must subtract half the\n"
                     "  dissociation energy.")
        if len(validos) > 1:
            segundo = validos[1][1]
            L.append(f"  Difference from the second site: "
                     f"{abs(segundo - mejor_e):.4f} eV")
            if abs(segundo - mejor_e) < 0.05:
                L.append("  The first two are within 50 meV: at this "
                         "precision one cannot\n  say which wins; finer "
                         "cutoffs and mesh are needed to separate them.")

    if run.slab_ok is False or run.mol_ok is False:
        L.append("")
        L.append("WARNING: a reference did not converge; E_ads inherits that error.")
    sin_conv = [run.sitios[i].etiqueta for i in range(len(run.sitios))
                if run.converged[i] is False]
    if sin_conv:
        L.append(f"NOT CONVERGED: {', '.join(sin_conv)}")
    return "\n".join(L)


def export(run: AdsorbRun, outdir: str = ".") -> list:
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    f = out / "ADSORCION.dat"
    eads = run.energias_ads
    lines = [provenance.header(
        f"adsorption energies of {run.molecula}",
        {"E_slab_eV": run.E_slab, "E_mol_eV": run.E_mol,
         "atomos_losa": run.natoms_slab}),
        f"# {'site':<12s} {'type':<8s} {'E_ads(eV)':>12s} "
        f"{'height(A)':>11s} {'contact(A)':>12s}"]
    nan = float("nan")
    for i, s in enumerate(run.sitios):
        if eads[i] is None:
            continue
        lines.append(
            f"  {s.etiqueta:<12s} {s.tipo:<8s} {eads[i]:>12.5f} "
            f"{(run.alturas[i] if run.alturas[i] is not None else nan):>11.3f} "
            f"{(run.contactos[i] if run.contactos[i] is not None else nan):>12.3f}")
    f.write_text("\n".join(lines) + "\n", encoding="utf-8")
    txt = out / "ADSORCION.txt"
    txt.write_text(report(run) + "\n", encoding="utf-8")
    return [str(f), str(txt)]


def plot(run: AdsorbRun, outfile: str = "adsorcion", formats="pdf,png",
         theme: str = None, size: str = None, family: str = None,
         background: str = None, palette=None, usetex: bool = None,
         width="single", journal: str = "generic", aspect: float = 0.70,
         mono: bool = False, dpi: int = None) -> list:
    """Barras de E_ads por sitio, ordenadas y coloreadas por tipo."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:                              # pragma: no cover
        raise RuntimeError("matplotlib is not installed.") from exc

    eads = run.energias_ads
    datos = [(run.sitios[i].etiqueta, run.sitios[i].tipo, eads[i])
             for i in range(len(run.sitios)) if eads[i] is not None]
    if not datos:
        raise FaltanDatos("there are no adsorption energies to plot.")
    datos.sort(key=lambda t: t[2])

    st = qstyle.apply(theme, size=size, family=family, background=background,
                      palette=palette, usetex=usetex, mono=mono)
    fig, ax = qstyle.new_figure(width, journal, aspect)
    cols = qstyle.palette(3, mono=mono)
    color_de = {"top": cols[0], "bridge": cols[1], "hollow": cols[2]}
    x = np.arange(len(datos))
    ax.bar(x, [d[2] for d in datos], width=0.68,
           color=[color_de.get(d[1], cols[0]) for d in datos])
    ax.axhline(0.0, color=qstyle.INK_FAINT, lw=st["axis_line"])
    ax.set_xticks(x)
    ax.set_xticklabels([d[0] for d in datos], rotation=45, ha="right",
                       fontsize=st["legend"])
    ax.set_ylabel(r"$E_{\mathrm{ads}}$ (eV)")
    vistos = []
    for t in ("top", "bridge", "hollow"):
        if any(d[1] == t for d in datos):
            vistos.append(plt.Rectangle((0, 0), 1, 1, color=color_de[t], label=t))
    if len(vistos) > 1:
        ax.legend(handles=vistos, frameon=False, fontsize=st["legend"])
    written = qstyle.save(fig, outfile, formats, dpi=dpi, modulo="adsorcion")
    plt.close(fig)
    return written
