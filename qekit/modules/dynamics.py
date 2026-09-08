# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Análisis de trayectorias de dinámica molecular.

De una corrida de `pw.x calculation='md'` (o de cp.x) sale una lista de
posiciones contra el tiempo. Sola no dice nada; lo que informa son tres
funciones que se extraen de ella:

- **g(r), la función de distribución radial.** Cuántos vecinos hay a cada
  distancia. En un líquido o un amorfo es LA medida de estructura, y es
  directamente comparable con lo que sale de difracción de neutrones o
  rayos X. El primer pico da la distancia de enlace; su área, el número
  de coordinación.

- **MSD, el desplazamiento cuadrático medio.** Si crece linealmente, hay
  difusión, y la pendiente da el coeficiente D por la relación de
  Einstein: MSD = 6Dt en tres dimensiones. Si se aplana, el átomo está
  vibrando alrededor de un sitio y NO difunde.

- **VDOS, la densidad de estados vibracional.** La transformada de
  Fourier de la autocorrelación de velocidades. Da el espectro
  vibracional sin DFPT y, a diferencia de los fonones armónicos, INCLUYE
  la anarmonicidad y la temperatura.

LO QUE HAY QUE MIRAR ANTES DE CREERSE NINGUNA DE LAS TRES
---------------------------------------------------------
1. **El equilibrado.** Los primeros pasos de una MD no son la trayectoria
   de equilibrio: el sistema está soltando la energía de la configuración
   inicial. Analizarlos mete un sesgo. Olla-DFT descarta un tramo inicial
   (`--skip`) y avisa si la temperatura todavía tiene tendencia.

2. **La longitud.** Un D sacado de 2 ps de trayectoria no es un D. La
   regla práctica es que el MSD tiene que ser lineal durante al menos un
   orden de magnitud en tiempo, y que el desplazamiento supere varias
   veces la distancia interatómica. Olla-DFT calcula el ajuste solo en el
   tramo donde es lineal y reporta el R².

3. **g(r) y la caja.** Con condiciones periódicas, g(r) solo tiene
   sentido hasta la mitad de la arista menor de la celda. Más allá, la
   normalización deja de valer. Olla-DFT corta ahí.
"""

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from qekit.core import provenance
from qekit.core import style as qstyle
from qekit.core.compat import trapezoid
from qekit.core.errors import ErrorDeUso

BOHR_A = 0.529177210903
KB_RY = 6.33362e-6          # k_B en Ry/K


@dataclass
class Trayectoria:
    simbolos: list = field(default_factory=list)
    posiciones: np.ndarray = None     # (nsteps, nat, 3) en Å
    celda: np.ndarray = None          # (3,3) Å  (se asume constante)
    tiempos: np.ndarray = None        # fs
    temperaturas: np.ndarray = None   # K
    energias: np.ndarray = None       # eV (energía total)
    dt: float = None                  # fs entre pasos guardados

    @property
    def nsteps(self) -> int:
        return 0 if self.posiciones is None else len(self.posiciones)

    @property
    def natoms(self) -> int:
        return len(self.simbolos)


@dataclass
class AnalisisMD:
    tray: Trayectoria = None
    r: np.ndarray = None
    gr: np.ndarray = None
    gr_pares: dict = field(default_factory=dict)      # ("Si","O") -> g(r)
    coordinacion: dict = field(default_factory=dict)  # par -> (r_min, N)
    t_msd: np.ndarray = None
    msd: np.ndarray = None
    msd_especie: dict = field(default_factory=dict)
    D: dict = field(default_factory=dict)             # especie -> cm^2/s
    r2: dict = field(default_factory=dict)
    frecuencias: np.ndarray = None    # cm^-1
    vdos: np.ndarray = None
    T_media: float = None
    T_desv: float = None
    equilibrado: int = 0
    avisos: list = field(default_factory=list)


# ----------------------------------------------------------------------
# Lectura de la salida de pw.x
# ----------------------------------------------------------------------
_RE_TEMP = re.compile(r"temperature\s*=\s*(-?[\d.]+)\s*K")
_RE_ETOT = re.compile(r"^!\s+total energy\s*=\s*(-?[\d.]+)\s*Ry", re.M)
# pw.x lo dice ya convertido: "Time step = 20.00 a.u.,  0.9676 femto-seconds".
# Leer el numero en unidades atomicas y convertirlo a mano introduce un
# error de conversion gratuito, y ademas la etiqueta "dt =" del input no
# aparece en la salida.
_RE_DT = re.compile(r"Time step\s*=\s*[\d.]+\s*a\.u\.,\s*"
                    r"([\d.]+)\s*femto-seconds", re.I)


def leer_md(path, skip: int = 0) -> Trayectoria:
    """Lee posiciones, temperatura y energía de una salida de pw.x md."""
    p = Path(path)
    if p.is_dir():
        cand = sorted(list(p.glob("*.out")) + list(p.glob("md.out")))
        if not cand:
            raise ErrorDeUso(
                f"there is no .out output in {p}. This module reads the "
                "text output of pw.x with calculation='md'.")
        p = cand[0]
    texto = p.read_text(errors="ignore")
    if "ATOMIC_POSITIONS" not in texto:
        raise ErrorDeUso(
            f"{p.name} contains no ATOMIC_POSITIONS block: it does not look like "
            "the output of a molecular dynamics run. A normal scf does not move the "
            "atoms.")

    celda = _leer_celda(texto)
    marcos, simbolos = _leer_marcos(texto, celda, alat_de(texto))
    if not marcos:
        raise ErrorDeUso(f"could not read any configuration from {p.name}.")

    temps = np.array([float(m) for m in _RE_TEMP.findall(texto)], dtype=float)
    energias = np.array([float(m) * 13.605693122994
                         for m in _RE_ETOT.findall(texto)], dtype=float)
    mdt = _RE_DT.search(texto)
    dt_fs = float(mdt.group(1)) if mdt else 1.0

    pos = np.array(marcos, dtype=float)
    n = len(pos)
    tray = Trayectoria(simbolos=simbolos, celda=celda, dt=dt_fs)
    tray.posiciones = pos[skip:]
    tray.tiempos = np.arange(len(tray.posiciones), dtype=float) * dt_fs
    tray.temperaturas = temps[skip:skip + len(tray.posiciones)] \
        if len(temps) >= n else temps[skip:]
    tray.energias = energias[skip:skip + len(tray.posiciones)] \
        if len(energias) >= n else energias[skip:]
    return tray


def alat_de(texto):
    """alat en angstrom, de la cabecera de pw.x."""
    m = re.search(r"lattice parameter \(alat\)\s*=\s*([\d.]+)", texto)
    return float(m.group(1)) * BOHR_A if m else None


def _leer_celda(texto):
    a0 = alat_de(texto)
    # OJO: el patron tiene que anclarse en "a(i) = (...)". Buscar cualquier
    # parentesis captura tambien el "(1)" de la propia etiqueta a(1) y
    # devuelve seis grupos en vez de tres.
    filas = re.findall(r"a\(\d\)\s*=\s*\(([^)]*)\)", texto)
    if len(filas) >= 3 and a0:
        v = np.array([[float(x) for x in f.split()] for f in filas[:3]])
        return v * a0
    m2 = re.search(r"CELL_PARAMETERS[^\n]*\n((?:[^\n]*\n){3})", texto)
    if m2:
        return np.array([[float(x) for x in ln.split()]
                         for ln in m2.group(1).strip().split("\n")])
    raise ErrorDeUso("could not read the cell from the pw.x output.")


def _a_angstrom(pos, unidad, celda, alat):
    """Lleva las posiciones a angstrom, sean cuales sean sus unidades."""
    if unidad in ("angstrom", "angstroms"):
        return pos
    if unidad == "bohr":
        return pos * BOHR_A
    if unidad == "alat":
        if alat is None:
            raise ErrorDeUso(
                "the positions are in alat units but the output does not "
                "contain the value of alat.")
        return pos * alat
    if unidad in ("crystal", "crystal_sg"):
        if celda is None:
            raise ErrorDeUso(
                "the positions are in crystal coordinates but the cell "
                "could not be read.")
        return pos @ np.asarray(celda, dtype=float)
    raise ErrorDeUso(
        "unrecognized position unit '%s'. Olla-DFT understands angstrom, "
        "bohr, alat and crystal." % unidad)


def _leer_marcos(texto, celda=None, alat=None):
    """Configuraciones de la trayectoria, SIEMPRE devueltas en angstrom.

    pw.x escribe las posiciones en las unidades que se le pidieron: alat,
    bohr, angstrom o crystal. Tratarlas todas como angstrom es el error
    que hace que un g(r) salga a la escala equivocada sin que nada se
    queje.
    """
    bloques = re.findall(
        r"ATOMIC_POSITIONS\s*[({]?\s*(\w+)\s*[)}]?\s*\n"
        r"((?:\s*\w+[^\n]*\n)+)", texto)
    marcos, simbolos, unidades = [], [], []
    for unidad, cuerpo in bloques:
        filas, syms = [], []
        for ln in cuerpo.strip().split("\n"):
            partes = ln.split()
            if len(partes) < 4:
                continue
            try:
                xyz = [float(partes[1]), float(partes[2]), float(partes[3])]
            except ValueError:
                continue
            syms.append(partes[0]); filas.append(xyz)
        if not filas:
            continue
        if not simbolos:
            simbolos = syms
        if len(filas) != len(simbolos):
            continue
        marcos.append(_a_angstrom(np.array(filas, dtype=float),
                                  unidad.lower(), celda, alat))
        unidades.append(unidad.lower())
    if unidades and len(set(unidades)) > 1:
        raise ErrorDeUso(
            "the trajectory mixes position units (%s); Olla-DFT does not know "
            "which one to apply to each frame." % ", ".join(sorted(set(unidades))))
    return marcos, simbolos



# ----------------------------------------------------------------------
# g(r)
# ----------------------------------------------------------------------
def rdf(tray: Trayectoria, rmax: float = None, nbins: int = 200,
        pares=None) -> tuple:
    """Función de distribución radial, total y por pares de especies.

    El corte por omisión es la mitad de la arista menor de la celda: más
    allá, la esfera de radio r ya no cabe entera en la caja y la
    normalización por el volumen del cascarón deja de ser válida.
    """
    if tray.nsteps == 0:
        raise ErrorDeUso("the trajectory is empty.")
    celda = np.asarray(tray.celda, dtype=float)
    inv = np.linalg.inv(celda)
    L = min(np.linalg.norm(celda[i]) for i in range(3))
    if rmax is None:
        rmax = 0.5 * L
    rmax = min(rmax, 0.5 * L)

    bordes = np.linspace(0.0, rmax, nbins + 1)
    r = 0.5 * (bordes[1:] + bordes[:-1])
    dr = bordes[1] - bordes[0]
    V = abs(np.linalg.det(celda))
    nat = tray.natoms
    simbolos = np.array(tray.simbolos)
    especies = sorted(set(tray.simbolos))

    total = np.zeros(nbins)
    por_par = {}
    if pares is None:
        pares = [(a, b) for i, a in enumerate(especies)
                 for b in especies[i:]]
    for par in pares:
        por_par[tuple(sorted(par))] = np.zeros(nbins)

    for marco in tray.posiciones:
        d = marco[:, None, :] - marco[None, :, :]
        f = d @ inv
        f -= np.round(f)
        d = f @ celda
        dist = np.linalg.norm(d, axis=-1)
        iu = np.triu_indices(nat, k=1)
        dd = dist[iu]
        total += np.histogram(dd, bins=bordes)[0]
        s1, s2 = simbolos[iu[0]], simbolos[iu[1]]
        for par in por_par:
            sel = (((s1 == par[0]) & (s2 == par[1])) |
                   ((s1 == par[1]) & (s2 == par[0])))
            por_par[par] += np.histogram(dd[sel], bins=bordes)[0]

    cascaron = 4.0 * np.pi * r ** 2 * dr
    norm_total = tray.nsteps * (nat * (nat - 1) / 2.0) / V * cascaron
    gr = np.divide(total, norm_total, out=np.zeros_like(total),
                   where=norm_total > 0)

    conteos = {e: int(np.sum(simbolos == e)) for e in especies}
    gr_pares = {}
    for par, h in por_par.items():
        na, nb = conteos[par[0]], conteos[par[1]]
        npares = na * nb if par[0] != par[1] else na * (na - 1) / 2.0
        if npares <= 0:
            continue
        nrm = tray.nsteps * npares / V * cascaron
        gr_pares[par] = np.divide(h, nrm, out=np.zeros_like(h),
                                  where=nrm > 0)
    return r, gr, gr_pares


def coordinacion(r, gr, tray: Trayectoria, par=None) -> tuple:
    """(r del primer mínimo, número de vecinos hasta ahí).

    El número de coordinación es la integral de 4*pi*r^2*rho*g(r) hasta
    el primer mínimo. Elegir ese mínimo a ojo es lo habitual y lo que
    hace que dos personas obtengan números distintos; aquí se toma el
    primer mínimo local después del primer máximo, y se reporta cuál fue.
    """
    gr = np.asarray(gr, dtype=float)
    if gr.max() <= 0:
        return float("nan"), float("nan")
    imax = int(np.argmax(gr))
    imin = None
    for i in range(imax + 1, len(gr) - 1):
        if gr[i] <= gr[i - 1] and gr[i] <= gr[i + 1] and gr[i] < 1.0:
            imin = i
            break
    if imin is None:
        return float("nan"), float("nan")
    V = abs(np.linalg.det(np.asarray(tray.celda, dtype=float)))
    rho = tray.natoms / V
    n = trapezoid(4.0 * np.pi * r[:imin + 1] ** 2 * rho * gr[:imin + 1],
                  r[:imin + 1])
    return float(r[imin]), float(n)


# ----------------------------------------------------------------------
# MSD y difusión
# ----------------------------------------------------------------------
def desdoblar(tray: Trayectoria) -> np.ndarray:
    """Quita los saltos de la imagen periódica.

    Sin esto, un átomo que cruza la frontera de la caja aparece saltando
    de un lado al otro y el MSD sale disparado. Es el error clásico de
    calcular difusión desde una salida de MD.
    """
    celda = np.asarray(tray.celda, dtype=float)
    inv = np.linalg.inv(celda)
    frac = tray.posiciones @ inv
    saltos = np.zeros_like(frac)
    saltos[1:] = np.cumsum(np.round(np.diff(frac, axis=0)), axis=0)
    return (frac - saltos) @ celda


def msd(tray: Trayectoria, especie: str = None) -> tuple:
    """Desplazamiento cuadrático medio contra el retardo, con media móvil
    sobre todos los orígenes de tiempo (mucho menos ruidoso)."""
    pos = desdoblar(tray)
    if especie:
        sel = np.array(tray.simbolos) == especie
        if not sel.any():
            raise ErrorDeUso(f"there are no '{especie}' atoms in the trajectory.")
        pos = pos[:, sel, :]
    n = len(pos)
    retardos = np.arange(1, max(n // 2, 2))
    out = np.zeros(len(retardos))
    for i, tau in enumerate(retardos):
        d = pos[tau:] - pos[:-tau]
        out[i] = float(np.mean(np.sum(d ** 2, axis=-1)))
    return retardos * tray.dt, out


def difusion(t_fs, msd_A2, frac_ini: float = 0.2,
             frac_fin: float = 0.8) -> tuple:
    """D por Einstein (MSD = 6Dt), ajustando solo el tramo central.

    El principio de la curva es balístico (no difusivo) y el final tiene
    poquísimos orígenes de tiempo promediados, así que es puro ruido.
    Ajustar la curva entera es el error más común aquí.

    Devuelve (D en cm^2/s, R^2 del ajuste).
    """
    t = np.asarray(t_fs, dtype=float)
    y = np.asarray(msd_A2, dtype=float)
    i0, i1 = int(frac_ini * len(t)), int(frac_fin * len(t))
    if i1 - i0 < 3:
        return float("nan"), float("nan")
    tt, yy = t[i0:i1], y[i0:i1]
    A = np.vstack([tt, np.ones_like(tt)]).T
    (m, b), res, *_ = np.linalg.lstsq(A, yy, rcond=None)
    ss_tot = float(np.sum((yy - yy.mean()) ** 2))
    ss_res = float(res[0]) if len(res) else float(
        np.sum((yy - (m * tt + b)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    # m está en Å²/fs -> cm²/s : (1e-8 cm)^2 / (1e-15 s) = 1e-1
    D = m / 6.0 * 1e-1
    return float(D), float(r2)


# ----------------------------------------------------------------------
# VDOS
# ----------------------------------------------------------------------
def vdos(tray: Trayectoria, ventana: bool = True) -> tuple:
    """Densidad de estados vibracional por autocorrelación de velocidades.

    Las velocidades salen de derivar las posiciones DESDOBLADAS: con las
    posiciones tal cual, un cruce de frontera mete un pico enorme que se
    reparte por todo el espectro.
    """
    if tray.nsteps < 8:
        raise ErrorDeUso(
            f"more steps are needed for a VDOS (there are {tray.nsteps}). With "
            "fewer than a few hundred the frequency resolution is useless.")
    pos = desdoblar(tray)
    v = np.gradient(pos, tray.dt, axis=0)          # Å/fs
    n = len(v)
    # autocorrelación por FFT, sumada sobre átomos y componentes
    nfft = 1
    while nfft < 2 * n:
        nfft *= 2
    acf = np.zeros(n)
    F = np.fft.rfft(v, n=nfft, axis=0)
    S = np.sum(np.abs(F) ** 2, axis=(1, 2))
    c = np.fft.irfft(S, n=nfft)[:n]
    acf = c / np.arange(n, 0, -1)
    if acf[0] > 0:
        acf = acf / acf[0]
    if ventana:
        acf = acf * np.hanning(2 * n)[n:]

    espectro = np.abs(np.fft.rfft(acf, n=4 * n))
    freq_fs = np.fft.rfftfreq(4 * n, d=tray.dt)    # 1/fs
    # 1/fs -> cm^-1 :  nu[cm^-1] = f[1/fs] * 1e15 / c[cm/s]
    cm1 = freq_fs * 1e15 / 2.99792458e10
    return cm1, espectro


# ----------------------------------------------------------------------
# Análisis completo
# ----------------------------------------------------------------------
def analizar(tray: Trayectoria, rmax: float = None, nbins: int = 200,
             equilibrado: int = 0) -> AnalisisMD:
    a = AnalisisMD(tray=tray, equilibrado=equilibrado)
    a.r, a.gr, a.gr_pares = rdf(tray, rmax=rmax, nbins=nbins)
    for par, g in a.gr_pares.items():
        a.coordinacion[par] = coordinacion(a.r, g, tray)
    a.t_msd, a.msd = msd(tray)
    D, r2 = difusion(a.t_msd, a.msd)
    a.D["total"] = D
    a.r2["total"] = r2
    for e in sorted(set(tray.simbolos)):
        t, m = msd(tray, especie=e)
        a.msd_especie[e] = m
        a.D[e], a.r2[e] = difusion(t, m)
    try:
        a.frecuencias, a.vdos = vdos(tray)
    except ErrorDeUso as exc:
        a.avisos.append(f"No VDOS: {exc}")

    if tray.temperaturas is not None and len(tray.temperaturas) > 2:
        T = np.asarray(tray.temperaturas, dtype=float)
        a.T_media = float(T.mean())
        a.T_desv = float(T.std())
        mitad = len(T) // 2
        deriva = float(T[mitad:].mean() - T[:mitad].mean())
        if abs(deriva) > 0.15 * max(a.T_media, 1.0):
            a.avisos.append(
                f"The temperature drifted {deriva:+.0f} K between the first and "
                f"the second half\nof the trajectory (mean {a.T_media:.0f} "
                "K). That means it is still equilibrating:\ndiscard more "
                "steps with --skip before analysing.")
    return a


def report(a: AnalisisMD) -> str:
    t = a.tray
    lines = ["--- Molecular dynamics ---",
             f"Trajectory: {t.nsteps} configurations, {t.natoms} atoms, "
             f"dt = {t.dt:.2f} fs",
             f"Analysed duration: {t.nsteps * t.dt / 1000:.2f} ps"
             + (f"  ({a.equilibrado} equilibration steps were discarded)"
                if a.equilibrado else "")]
    if a.T_media is not None:
        lines.append(f"Temperature: {a.T_media:.0f} ± {a.T_desv:.0f} K")
    if t.nsteps * t.dt < 2000:
        lines.append(
            "  WARNING: less than 2 ps. Enough to see the structure, not for "
            "a diffusion\n  coefficient: at least one order of "
            "magnitude more would be needed.")

    lines += ["", "Structure — g(r):"]
    rmax = float(a.r[-1]) if a.r is not None else 0.0
    if a.gr is not None and float(np.max(a.gr)) <= 0:
        lines += [f"  g(r) is empty up to the cutoff of {rmax:.2f} Å.",
                  "  That cutoff is half the cell edge, which is as far "
                  "as the normalization\n  of g(r) makes sense with "
                  "periodic boundary conditions. If the first bond\n  "
                  "distance lies above it, the cell is too small "
                  "to extract\n  structure from it: build a supercell."]
    elif a.coordinacion:
        lines.append(f"  {'pair':>10s} {'1st peak':>10s} {'1st min':>9s} "
                     f"{'coordination':>13s}")
        for par, (rmin, ncoord) in sorted(a.coordinacion.items()):
            g = a.gr_pares[par]
            rpico = a.r[int(np.argmax(g))] if g.max() > 0 else float("nan")
            lines.append(f"  {'-'.join(par):>10s} {rpico:9.3f} Å "
                         f"{rmin:8.3f} Å {ncoord:12.2f}")
        lines.append("  The coordination is integrated up to the FIRST MINIMUM of "
                     "g(r); that cutoff is\n  a convention, not a measurement.")

    lines += ["", "Diffusion — MSD versus lag time:"]
    lines.append(f"  {'species':>10s} {'D (cm²/s)':>14s} {'R²':>8s}")
    for e, D in sorted(a.D.items()):
        marca = ""
        if not np.isnan(a.r2.get(e, np.nan)) and a.r2[e] < 0.95:
            marca = "   <- the MSD is NOT linear: no diffusion, or not enough time"
        lines.append(f"  {e:>10s} {D:14.3e} {a.r2.get(e, float('nan')):8.4f}"
                     f"{marca}")
    lines.append("  D comes from MSD = 6Dt fitting only the central segment: the "
                 "beginning is\n  ballistic and the end has too little "
                 "averaging.")

    if a.frecuencias is not None and a.vdos is not None:
        v = a.vdos / max(a.vdos.max(), 1e-30)
        picos = [(a.frecuencias[i], v[i]) for i in range(1, len(v) - 1)
                 if v[i] > v[i - 1] and v[i] >= v[i + 1] and v[i] > 0.15]
        picos.sort(key=lambda p: -p[1])
        lines += ["", "Vibrations — VDOS (includes anharmonicity and "
                  "temperature):"]
        if picos:
            for f, alto in picos[:5]:
                lines.append(f"  {f:8.1f} cm⁻¹   intensity {alto:.2f}")
        lines.append(f"  Frequency resolution: "
                     f"{1e15 / (t.nsteps * t.dt * 2.99792458e10):.1f} cm⁻¹ "
                     "(set by the duration of the trajectory).")

    for w in a.avisos:
        lines += ["", w]
    return "\n".join(lines)


def export(a: AnalisisMD, outdir: str = ".") -> list:
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    escritos = []
    cab = provenance.header_plain(
        "molecular dynamics",
        {"pasos": a.tray.nsteps, "dt_fs": a.tray.dt,
         "T_media_K": None if a.T_media is None else round(a.T_media, 1)},
        titulo="Trajectory analysis")

    f = out / "MD_RDF.dat"
    cols = [a.r, a.gr] + [a.gr_pares[p] for p in sorted(a.gr_pares)]
    nombres = ["r(A)", "g_total"] + [f"g_{'-'.join(p)}"
                                     for p in sorted(a.gr_pares)]
    np.savetxt(f, np.column_stack(cols), fmt="%12.6f",
               header=cab + "\n" + "  ".join(f"{n:>12s}" for n in nombres),
               comments="# ")
    escritos.append(str(f))

    f = out / "MD_MSD.dat"
    cols = [a.t_msd, a.msd] + [a.msd_especie[e] for e in sorted(a.msd_especie)]
    nombres = ["t(fs)", "MSD(A2)"] + [f"MSD_{e}"
                                      for e in sorted(a.msd_especie)]
    np.savetxt(f, np.column_stack(cols), fmt="%14.6f",
               header=cab + "\n" + "  ".join(f"{n:>14s}" for n in nombres),
               comments="# ")
    escritos.append(str(f))

    if a.vdos is not None:
        f = out / "MD_VDOS.dat"
        np.savetxt(f, np.column_stack([a.frecuencias, a.vdos]), fmt="%14.6f",
                   header=cab + "\n   frequency(cm-1)           VDOS",
                   comments="# ")
        escritos.append(str(f))

    txt = out / "MD.txt"
    txt.write_text(report(a) + "\n", encoding="utf-8")
    escritos.append(str(txt))
    return escritos


def plot(a: AnalisisMD, outfile: str = "md", formats="pdf,png",
         theme: str = None, family: str = None, background: str = None,
         palette=None, usetex: bool = None, width="double",
         journal: str = "generic", mono: bool = False,
         dpi: int = None) -> list:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:                          # pragma: no cover
        raise RuntimeError("matplotlib is not installed.") from exc

    st = qstyle.apply(theme, family=family, background=background,
                      palette=palette, usetex=usetex, mono=mono)
    ancho = qstyle.figure_size(width, journal, 0.32)
    fig, axes = plt.subplots(1, 3, figsize=ancho, layout="constrained")
    colores = qstyle.palette(6, mono=mono)

    ax = axes[0]
    for i, par in enumerate(sorted(a.gr_pares)):
        ax.plot(a.r, a.gr_pares[par], lw=1.1,
                label="–".join(par), **qstyle.style_line(i, 6, mono=mono))
    ax.axhline(1.0, color=qstyle.INK_FAINT, lw=st["axis_line"],
               dashes=[3.5, 2.0])
    ax.set_xlabel("r (Å)"); ax.set_ylabel("g(r)")
    ax.set_xlim(0, a.r[-1]); ax.set_ylim(bottom=0)
    ax.legend(frameon=False, fontsize="small")

    ax = axes[1]
    for i, e in enumerate(sorted(a.msd_especie)):
        ax.plot(a.t_msd / 1000.0, a.msd_especie[e], lw=1.1, label=e,
                **qstyle.style_line(i, 6, mono=mono))
    ax.set_xlabel("lag time (ps)"); ax.set_ylabel(r"MSD (Å$^2$)")
    ax.set_xlim(left=0); ax.set_ylim(bottom=0)
    ax.legend(frameon=False, fontsize="small")

    ax = axes[2]
    if a.vdos is not None:
        sel = a.frecuencias < 2000
        ax.plot(a.frecuencias[sel], a.vdos[sel] / a.vdos[sel].max(), lw=1.1,
                color=colores[0])
    ax.set_xlabel(r"$\tilde\nu$ (cm$^{-1}$)"); ax.set_ylabel("VDOS (norm.)")
    ax.set_ylim(bottom=0)
    for i, ax in enumerate(axes):
        qstyle.panel_label(ax, f"({'abc'[i]})")
    return qstyle.save(fig, outfile, formats, dpi=dpi)
