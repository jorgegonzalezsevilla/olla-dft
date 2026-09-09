# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Validación cruzada: la misma cantidad por caminos independientes.

Cada módulo de Olla-DFT se valida contra la literatura, pero eso no detecta
un error que afecte a un módulo entero de forma sistemática. Lo que sí lo
detecta es calcular la MISMA cantidad por dos rutas físicamente
independientes y compararlas:

    modulo volumetrico   ->  ajuste de la EOS   vs   traza de las Cij
    velocidad del sonido ->  Cij                vs   pendiente acustica
    temperatura de Debye ->  velocidades        vs   DOS de fonones
    gap                  ->  estructura de bandas vs extrapolacion de Tauc
    C_v a T alta         ->  DOS de fonones     vs   limite de Dulong-Petit
    numero de modos      ->  integral de la DOS vs   3N

No cuesta ningun calculo nuevo: son resultados que ya estan en disco.

SOBRE LAS TOLERANCIAS
---------------------
Cada cruce lleva la suya, y no son arbitrarias:

- B0 por dos rutas es la MISMA cantidad y debe coincidir al ~5 %;
- las velocidades del sonido tambien, pero la pendiente acustica en q->0
  es justo lo que peor interpola una malla de q gruesa, asi que un
  desacuerdo ahi acusa a la malla antes que al modulo elastico;
- las dos temperaturas de Debye NO son la misma definicion (una es el
  limite acustico, la otra usa todo el espectro): ahi la tolerancia es
  amplia a proposito, y coincidir al 1 % seria sospechoso, no bueno.

Un cruce que falla no dice cual de los dos caminos esta mal. Por eso cada
uno lleva un diagnostico de que mirar primero.
"""

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

KB_EV = 8.617333262e-5


@dataclass
class Check:
    nombre: str = ""
    ruta_a: str = ""
    valor_a: float = None
    ruta_b: str = ""
    valor_b: float = None
    unidad: str = ""
    tolerancia: float = 0.05
    diagnostico: str = ""

    @property
    def desvio(self):
        """Desviación RELATIVA, salvo que la referencia sea cero.

        Hay cruces cuya respuesta correcta ES cero —la fase de Berry de un
        cristal centrosimétrico, por ejemplo— y ahí no existe desviación
        relativa. Antes se devolvía None, el cruce se daba por "sin datos" y
        el informe reventaba al formatearlo. Con referencia cero la
        desviación es el valor absoluto, y la tolerancia se lee en las
        unidades del cruce.
        """
        if self.valor_a is None or self.valor_b is None:
            return None
        if self.valor_a == 0:
            return abs(self.valor_b)
        return abs(self.valor_b - self.valor_a) / abs(self.valor_a)

    @property
    def relativa(self):
        return self.valor_a not in (None, 0)

    @property
    def ok(self):
        d = self.desvio
        return None if d is None else bool(d <= self.tolerancia)


@dataclass
class CrossResult:
    checks: list = field(default_factory=list)
    disponibles: list = field(default_factory=list)
    faltantes: list = field(default_factory=list)


# ----------------------------------------------------------------------
def _cargar(project: Path) -> dict:
    """Busca en la carpeta los resultados que Olla-DFT ya haya escrito."""
    p = Path(project)
    d = {}
    for nombre, patron in (("elastic", "ELASTIC_C.dat"),
                           ("eos", "EOS.txt"),
                           ("fonones_dos", "FONONES_DOS.dat"),
                           ("fonones_bandas", "FONONES_BANDAS.dat"),
                           ("optics", "OPTICS.dat"),
                           ("kappa", "KAPPA.dat"),
                           ("berry", "BERRY.dat"),
                           ("wannier", "WANNIER_centros.dat"),
                           ("esm", "ESM.dat"),
                           ("wf", "WF.dat"),
                           ("strain", "STRAIN.dat")):
        hits = list(p.rglob(patron))
        if hits:
            d[nombre] = hits[0]
    return d


def _leer_eos_b0(path) -> float:
    for linea in Path(path).read_text(errors="ignore").splitlines():
        if "B0" in linea and "GPa" in linea:
            for tok in linea.replace("=", " ").split():
                try:
                    v = float(tok)
                    if 0.1 < v < 1e4:
                        return v
                except ValueError:
                    continue
    return None


def _leer_cabecera(path, clave):
    """Un valor escrito como '# clave = valor' en la cabecera de un .dat."""
    for linea in Path(path).read_text(errors="ignore").splitlines():
        if not linea.startswith("#") or clave not in linea:
            continue
        try:
            return float(linea.split("=")[1].split()[0])
        except (IndexError, ValueError):
            continue
    return None


def _leer_kappa_300(path):
    """κ medio a la temperatura más cercana a 300 K, en W/m/K."""
    d = np.loadtxt(path, comments="#")
    if d.ndim == 1:
        d = d[None, :]
    i = int(np.argmin(np.abs(d[:, 0] - 300.0)))
    return float(d[i, -1]), float(d[i, 0])


def _leer_berry_el(path):
    """Fase electrónica de Berry del punto de carga cero (unidades de QE)."""
    d = np.loadtxt(path, comments="#")
    if d.ndim == 1:
        d = d[None, :]
    i = int(np.argmin(np.abs(d[:, 0])))          # lambda o carga más cercana a 0
    return float(d[i, 2])


def _fase_de_centros(path, cell, gdir=3, spin=2.0):
    """La misma fase, desde los centros de Wannier: −f·Σ_n (r̄_n·b)/2π."""
    d = np.loadtxt(path, comments="#")
    if d.ndim == 1:
        d = d[None, :]
    frac = d[:, 1:4] @ np.linalg.inv(np.asarray(cell, float))
    return float(-spin * frac[:, int(gdir) - 1].sum())


def _leer_cij(path) -> np.ndarray:
    datos = np.loadtxt(path, comments="#")
    return datos if datos.shape == (6, 6) else None


# ----------------------------------------------------------------------
def run(project=".", masas=None, volumen=None, natoms=None,
        gap_bandas: float = None, C: np.ndarray = None,
        b0_eos: float = None, qdist=None, band_freqs=None,
        dos_w=None, dos=None, gap_tauc: float = None,
        cell=None, n_primitiva: int = None) -> CrossResult:
    """Ejecuta todos los cruces para los que haya datos."""
    from qekit.modules import derived

    res = CrossResult()
    encontrados = _cargar(project) if project else {}

    if C is None and "elastic" in encontrados:
        C = _leer_cij(encontrados["elastic"])
    if b0_eos is None and "eos" in encontrados:
        b0_eos = _leer_eos_b0(encontrados["eos"])
    if dos_w is None and "fonones_dos" in encontrados:
        datos = np.loadtxt(encontrados["fonones_dos"], comments="#")
        dos_w, dos = datos[:, 0], datos[:, 1]
    if qdist is None and "fonones_bandas" in encontrados:
        datos = np.loadtxt(encontrados["fonones_bandas"], comments="#")
        qdist, band_freqs = datos[:, 0], datos[:, 1:]

    for k, v in (("elastic constants", C is not None),
                 ("lattice κ (fc3)", "kappa" in encontrados),
                 ("Berry phase", "berry" in encontrados),
                 ("Wannier centres", "wannier" in encontrados),
                 ("work function (ESM)", "esm" in encontrados),
                 ("work function (planar potential)", "wf" in encontrados),
                 ("strain sweep", "strain" in encontrados),
                 ("equation of state", b0_eos is not None),
                 ("phonon DOS", dos_w is not None),
                 ("phonon dispersion", qdist is not None),
                 ("band gap", gap_bandas is not None),
                 ("Tauc gap", gap_tauc is not None)):
        (res.disponibles if v else res.faltantes).append(k)

    # --- 1. modulo volumetrico: EOS contra Cij -----------------------
    if C is not None and b0_eos is not None:
        from qekit.modules import elastic
        m = elastic.moduli(C)
        res.checks.append(Check(
            nombre="bulk modulus B₀",
            ruta_a="equation of state fit", valor_a=b0_eos,
            ruta_b="trace of the elastic constants",
            valor_b=m.B_hill, unidad="GPa", tolerancia=0.05,
            diagnostico=(
                "They are the SAME quantity by two routes. If they differ: check "
                "that the cell of the\nelastic constants was relaxed (low residual "
                "stress) and that the EOS has\nenough points on both sides of "
                "the minimum.")))

    # --- 2. velocidades del sonido -----------------------------------
    if (C is not None and qdist is not None and masas is not None
            and volumen):
        rho = derived.density(masas, volumen)
        dirs = derived.cubic_directional(C, rho)
        ac = derived.acoustic_velocities(qdist, band_freqs)
        if dirs and ac:
            res.checks.append(Check(
                nombre="longitudinal velocity [100]",
                ruta_a="elastic constants: √(C₁₁/ρ)",
                valor_a=dirs["v_l_100"],
                ruta_b="slope of the LA branch at Γ",
                valor_b=ac["v_l"], unidad="m/s", tolerancia=0.10,
                diagnostico=(
                    "The acoustic slope at q→0 is what a coarse q-grid "
                    "interpolates worst.\nIf it fails, suspect the grid "
                    "before the Cij.")))
            res.checks.append(Check(
                nombre="transverse velocity [100]",
                ruta_a="elastic constants: √(C₄₄/ρ)",
                valor_a=dirs["v_t_100"],
                ruta_b="slope of the TA branch at Γ",
                valor_b=ac["v_t1"], unidad="m/s", tolerancia=0.10,
                diagnostico=(
                    "The TRANSVERSE branches are the flattest and the ones "
                    "that come out worst from a small\nq-grid — in silicon "
                    "with 2x2x2 the error exceeds 40 %. If the\n"
                    "longitudinal one agrees and this one does not, it is the q-grid, not the "
                    "elastic constants.\nDensify it (4x4x4 or more) before trusting "
                    "either of the two.")))

    # --- 3. temperatura de Debye -------------------------------------
    if (C is not None and dos_w is not None and masas is not None
            and volumen and natoms):
        from qekit.modules import elastic
        m = elastic.moduli(C)
        rho = derived.density(masas, volumen)
        _vl, _vt, vm = derived.sound_velocities(m.B_hill, m.G_hill, rho)
        td_el = derived.debye_from_velocity(vm, natoms, volumen)
        td_dos = derived.debye_from_dos(dos_w, dos, natoms)
        if td_el and td_dos:
            res.checks.append(Check(
                nombre="Debye temperature",
                ruta_a="sound velocities (acoustic limit)",
                valor_a=td_el,
                ruta_b="second moment of the phonon DOS",
                valor_b=td_dos, unidad="K", tolerancia=0.30,
                diagnostico=(
                    "NOTE: they are NOT the same definition. The elastic one is the "
                    "low-temperature limit\n(acoustic only); the DOS one "
                    "uses the whole spectrum, optical modes included, and comes out\n"
                    "higher. They are cross-checked to catch a blunder, not to "
                    "coincide:\nagreement to 1 % would be suspicious.")))

    # --- 4. gap: bandas contra Tauc ----------------------------------
    if gap_bandas is not None and gap_tauc is not None:
        res.checks.append(Check(
            nombre="optical gap",
            ruta_a="band structure (direct gap)", valor_a=gap_bandas,
            ruta_b="Tauc extrapolation on α(E)", valor_b=gap_tauc,
            unidad="eV", tolerancia=0.06,
            diagnostico=(
                "epsilon.x does not include phonon-assisted transitions, "
                "so the absorption\nedge is the DIRECT gap, not the "
                "fundamental one. If you compare against the\nfundamental gap of an "
                "indirect semiconductor, the difference is physics, not an "
                "error.")))

    # --- 5. C_v a alta T contra Dulong-Petit -------------------------
    if dos_w is not None and natoms:
        cv = _cv_alta_T(dos_w, dos, natoms, T=1500.0)
        dp = 3.0 * natoms * KB_EV * 1000.0        # meV/K por celda
        if cv:
            res.checks.append(Check(
                nombre="C_v in the classical limit",
                ruta_a="Dulong–Petit: 3N·k_B", valor_a=dp,
                ruta_b="integral of the phonon DOS at 1500 K",
                valor_b=cv, unidad="meV/K per cell", tolerancia=0.03,
                diagnostico=(
                    "At high temperature every harmonic C_v tends to 3N·k_B. "
                    "If it falls short, the DOS\nis badly normalized or is missing "
                    "spectrum; if it overshoots, there are extra modes.")))

    # --- 6. numero de modos ------------------------------------------
    if dos_w is not None and natoms:
        from qekit.core.compat import trapezoid
        total = float(trapezoid(np.asarray(dos), np.asarray(dos_w)))
        res.checks.append(Check(
            nombre="number of modes",
            ruta_a="3N by construction", valor_a=3.0 * natoms,
            ruta_b="integral of the phonon DOS", valor_b=total,
            unidad="modes", tolerancia=0.05,
            diagnostico=(
                "The integral of the DOS must give exactly 3N. If not, "
                "the matdyn\ninterpolation grid is too poor or "
                "the frequency range cuts off\nspectrum.")))
    # --- 7. kappa de red: tercer orden contra el modelo de Slack ------
    if "kappa" in encontrados and C is not None and masas is not None \
            and volumen and natoms:
        from qekit.modules import elastic
        k300, T_usada = _leer_kappa_300(encontrados["kappa"])
        m = elastic.moduli(C)
        rho = derived.density(masas, volumen)
        _vl, _vt, vm = derived.sound_velocities(m.B_hill, m.G_hill, rho)
        td = derived.debye_from_velocity(vm, natoms, volumen)
        gam = derived.gruneisen_from_poisson(m.nu)
        # la n de Slack es la de la celda primitiva, no la de la celda dada
        ks = derived.slack(td, gam, float(np.mean(masas)), natoms, volumen,
                           T=T_usada,
                           n_celda=n_primitiva or natoms) if (td and gam) else None
        if ks:
            res.checks.append(Check(
                nombre="lattice thermal conductivity",
                ruta_a=f"phonon Boltzmann equation with fc3 "
                       f"({T_usada:.0f} K)",
                valor_a=k300,
                ruta_b="Slack model from the elastic constants",
                valor_b=ks, unidad="W/m/K", tolerancia=0.60,
                diagnostico=(
                    "The tolerance is 60 % ON PURPOSE: Slack is an "
                    "order-of-magnitude estimate\nwith an empirical prefactor, "
                    "not a calculation. It serves to detect that the\n"
                    "fc3 lacks convergence or that the sign of something is "
                    "wrong, not to fine-tune.\nIf they differ by a factor of 3, "
                    "suspect the fc3 supercell first.")))

    # --- 8. fase de Berry: lberry contra los centros de Wannier -------
    if "berry" in encontrados and "wannier" in encontrados and cell is not None:
        try:
            fa = _leer_berry_el(encontrados["berry"])
            fb = _fase_de_centros(encontrados["wannier"], cell)
            # las dos están definidas módulo 2: se comparan en la misma rama
            fb = fb - 2.0 * np.round((fb - fa) / 2.0)
            res.checks.append(Check(
                nombre="electronic Berry phase",
                ruta_a="lberry: determinant of overlaps along k-strings",
                valor_a=fa,
                ruta_b="Wannier centres: −2·Σ (r̄·b)/2π",
                valor_b=fb, unidad="(quantum = 2)", tolerancia=0.05,
                diagnostico=(
                    "They are the SAME Berry phase from two routines that do not "
                    "share a line of\ncode. Their agreement is the "
                    "strongest validation available here. If they disagree,\n"
                    "the first thing to check is that both use the same "
                    "direction (gdir) and the same\nk-point grid.")))
        except Exception:                                   # noqa: BLE001
            pass

    # --- 9. funcion trabajo: ESM contra la meseta del potencial -------
    if "esm" in encontrados and "wf" in encontrados:
        try:
            d = np.loadtxt(encontrados["esm"], comments="#")
            if d.ndim == 1:
                d = d[None, :]
            i = int(np.argmin(np.abs(d[:, 0])))     # la losa neutra
            phi_esm = float(d[i, 4])
            phi_wf = _leer_cabecera(encontrados["wf"], "Phi_eV")
            if phi_wf is not None:
                res.checks.append(Check(
                    nombre="work function",
                    ruta_a="ESM: the vacuum level is zero by "
                           "construction",
                    valor_a=phi_esm,
                    ruta_b="plateau of the planar potential from the pp.x cube",
                    valor_b=phi_wf, unidad="eV", tolerancia=0.05,
                    diagnostico=(
                        "With bc1 the ESM vacuum level is exactly zero and "
                        "there is no plateau to\nfit; the cube route does "
                        "fit it, and that is why it needs more vacuum. If\n"
                        "they differ, look at the flatness reported by the second: "
                        "it is almost always that the\nperiodic calculation "
                        "lacked vacuum, not that ESM is wrong.")))
        except Exception:                                   # noqa: BLE001
            pass

    # --- 10. modulo volumetrico: EOS contra la presion del barrido ----
    if "strain" in encontrados and b0_eos is not None:
        try:
            d = np.loadtxt(encontrados["strain"], comments="#")
            eps, P = d[:, 0], d[:, 3]
            bien = np.isfinite(P)
            if bien.sum() >= 3:
                # hidrostática: V = V0(1+ε)³  ->  B = −dP/d(lnV) = −dP/dε / 3
                pend = np.polyfit(eps[bien], P[bien], 1)[0]
                # La columna 3 de STRAIN.dat ya viene en GPa (la escribe
                # strain.export con cabecera P(GPa), desde res.pressure, que
                # qeout convierte con HA_BOHR3_GPA). No hay nada que pasar
                # de kbar: multiplicar por 0.1 dejaba B0 diez veces pequeño
                # y esta tercera ruta discrepaba siempre.
                b0_strain = -pend / 3.0
                if b0_strain > 0:
                    res.checks.append(Check(
                        nombre="bulk modulus B₀ (third route)",
                        ruta_a="equation of state fit",
                        valor_a=b0_eos,
                        ruta_b="slope of the pressure in the strain "
                               "sweep",
                        valor_b=b0_strain, unidad="GPa", tolerancia=0.10,
                        diagnostico=(
                            "Only valid if the sweep was HYDROSTATIC: with "
                            "biaxial or\nuniaxial strain the relation "
                            "between pressure and ε is different and this cross-check compares "
                            "apples\nwith oranges. Check it before trusting it.")))
        except Exception:                                   # noqa: BLE001
            pass

    return res


def _cv_alta_T(w, g, natoms, T=1500.0):
    """C_v armónica a temperatura alta, en meV/K por celda."""
    from qekit.core.compat import trapezoid
    w = np.asarray(w, dtype=float)
    g = np.asarray(g, dtype=float)
    m = w > 1.0
    w, g = w[m], g[m]
    if w.size < 3:
        return None
    norm = trapezoid(g, w)
    if norm <= 0:
        return None
    g = g * (3.0 * natoms / norm)
    e = w * 1.239841984e-4                      # cm^-1 -> eV
    x = e / (KB_EV * T)
    x = np.clip(x, 1e-9, 300.0)
    occ = 1.0 / np.expm1(x)
    cv = KB_EV * trapezoid(x ** 2 * np.exp(x) * occ ** 2 * g, w)
    return float(cv * 1000.0)


def report(res: CrossResult) -> str:
    lines = ["--- Cross-validation ---"]
    if res.disponibles:
        lines.append("Results found: " + ", ".join(res.disponibles))
    if res.faltantes:
        lines.append("Not available: " + ", ".join(res.faltantes))
    lines.append("")
    if not res.checks:
        lines.append(
            "There are no two independent routes to cross-check yet. Each cross-check "
            "needs TWO\nmodules: for example elastic constants + EOS, or "
            "elastic constants + phonons.")
        return "\n".join(lines)

    fallos = [c for c in res.checks if c.ok is False]
    lines.append(f"{len(res.checks)} cross-checks  |  "
                 f"{len(res.checks) - len(fallos)} agree  |  "
                 f"{len(fallos)} do NOT")
    lines.append("")
    for c in res.checks:
        marca = "OK  " if c.ok else ("FAIL " if c.ok is False else "  ?  ")
        if c.desvio is None:
            lines.append(f"[{marca}] {c.nombre}  (insufficient data)")
        elif c.relativa:
            lines.append(f"[{marca}] {c.nombre}  ({c.desvio * 100:.1f} % "
                         f"deviation, tolerance {c.tolerancia * 100:.0f} %)")
        else:
            lines.append(f"[{marca}] {c.nombre}  (must be zero; gives "
                         f"{c.desvio:.2e}, tolerance {c.tolerancia:g})")
        va = "—" if c.valor_a is None else f"{c.valor_a:.4g}"
        vb = "—" if c.valor_b is None else f"{c.valor_b:.4g}"
        lines.append(f"         {c.ruta_a}: {va} {c.unidad}")
        lines.append(f"         {c.ruta_b}: {vb} {c.unidad}")
        if c.ok is False:
            for l in c.diagnostico.splitlines():
                lines.append(f"         > {l}")
        lines.append("")

    if fallos:
        lines.append(
            "A failing cross-check does NOT say which of the two routes is wrong: "
            "it says that one of\nthe two is. The diagnostic of each one "
            "indicates what to look at first.")
    else:
        lines.append(
            "All cross-checks agree. This is the strongest evidence one "
            "can have without\nleaving the calculation itself: two "
            "independent routes do not make the same mistake by\nchance.")
    return "\n".join(lines)
