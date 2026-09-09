# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Transporte balístico: conductancia con `pwcond.x`.

QUÉ ES Y EN QUÉ SE DIFERENCIA DE `olla-dft transport`
--------------------------------------------------
`olla-dft transport` calcula transporte DIFUSIVO: electrones que se dispersan
muchas veces mientras cruzan el material, descritos por la ecuación de
Boltzmann. Es lo que corresponde a un cristal macroscópico.

Esto es lo contrario. En un nanocontacto, una molécula entre dos
electrodos o un nanohilo corto, el electrón cruza SIN dispersarse. Ahí no
hay conductividad: hay CONDUCTANCIA, y la da la fórmula de Landauer:

    G = G0 * T(E_F)        con G0 = 2e^2/h = 7.748e-5 S = 1/(12.906 kOhm)

T es la probabilidad de que un electrón que entra por la izquierda salga
por la derecha, sumada sobre todos los canales abiertos. Como T <= número
de canales, la conductancia viene CUANTIZADA en escalones de G0 — y ver
esos escalones es la comprobación de que el cálculo está bien.

LA GEOMETRÍA IMPORTA MÁS QUE NADA
---------------------------------
`pwcond.x` no acepta cualquier estructura. Necesita:

- el transporte a lo largo de **z**, siempre;
- un electrodo IZQUIERDO periódico en z, con su propio cálculo scf;
- opcionalmente una región de dispersión (la molécula, el defecto);
- un electrodo DERECHO, que en Olla-DFT es siempre el MISMO que el
  izquierdo (`ikind=1`). pwcond.x admite electrodos distintos (`ikind=2`,
  con `prefixr` y `bdr`), pero Olla-DFT no lo prepara: haría falta un
  tercer scf y la comprobación de que las tres celdas empalman.

Las regiones tienen que tener la MISMA celda en el plano xy, y los límites
(`bdl`, `bds`) se dan en unidades de alat a lo largo de z. Para pwcond.x el
electrodo ocupa de z = 0 a z = bdl y la región de dispersión de z = 0 a
z = bds, cada una en SU celda: bdl y bds son las LONGITUDES de esas celdas a
lo largo de z (celldm(3) si alat es |a1|), no la altura del último átomo.
Poner mal esos límites es la causa número uno de resultados sin sentido, y
Olla-DFT los calcula de la geometría en vez de dejarlos al ojo.

LOS DOS MODOS
-------------
- `ikind=0`: solo la **estructura de bandas compleja** del electrodo. Da
  el número de canales abiertos a cada energía, que es la COTA SUPERIOR de
  la conductancia. Es barato y es lo primero que hay que mirar.
- `ikind=1`: la conductancia de verdad, con la región de dispersión en
  medio y el mismo electrodo a los dos lados. Mucho más caro.
"""

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from qekit.core import provenance
from qekit.core import style as qstyle
from qekit.core.errors import ErrorDeUso
from qekit.modules import sweep

#: Cuanto de conductancia 2e^2/h, en siemens.
G0 = 7.748091729e-5
#: Su inverso, en ohmios.
R0 = 1.0 / G0


@dataclass
class CondRun:
    energias: np.ndarray = None      # eV respecto de E_F
    transmision: np.ndarray = None   # T(E), adimensional
    canales: np.ndarray = None       # nº de canales abiertos por energía
    ikind: int = None
    G_fermi: float = None            # en unidades de G0
    avisos: list = field(default_factory=list)

    @property
    def G_siemens(self) -> float:
        return (self.G_fermi or float("nan")) * G0

    @property
    def R_ohm(self) -> float:
        g = self.G_fermi
        return R0 / g if g else float("inf")


# ----------------------------------------------------------------------
# Geometría
# ----------------------------------------------------------------------
def limites_z(atoms, margen: float = 0.0) -> tuple:
    """(z mínimo, z máximo) de los átomos en unidades de alat.

    pwcond.x pide los límites de cada región en unidades de alat a lo
    largo de z, y alat es |a1|. Calcularlos de la geometría evita el error
    más común del módulo: cortar la región por donde no toca.
    """
    celda = np.array(atoms.get_cell(), dtype=float)
    alat = float(np.linalg.norm(celda[0]))
    z = atoms.get_positions()[:, 2]
    return (float(z.min() - margen) / alat, float(z.max() + margen) / alat)


def longitud_z(atoms) -> float:
    """Longitud de la celda a lo largo de z, en unidades de alat (= |a1|).

    Es lo que pwcond.x espera en `bdl` (electrodo) y `bds` (región de
    dispersión): la frontera derecha de cada región es el final de SU celda,
    no la altura del último átomo. Con un átomo en z = 0 y otro a mitad de
    celda, la altura máxima atómica sería medio periodo y pwcond.x
    empalmaría las regiones donde no toca.
    """
    celda = np.array(atoms.get_cell(), dtype=float)
    alat = float(np.linalg.norm(celda[0]))
    return float(np.linalg.norm(celda[2])) / alat


def comprobar_geometria(atoms) -> list:
    """Lo que tiene que cumplir una estructura para pwcond.x."""
    problemas = []
    celda = np.array(atoms.get_cell(), dtype=float)
    # el eje de transporte es z y tiene que ser ortogonal al plano
    if abs(celda[2, 0]) > 1e-6 or abs(celda[2, 1]) > 1e-6:
        problemas.append(
            "the third lattice vector is not parallel to z. pwcond.x ALWAYS "
            "transports along\nz, and the cell has to be oriented "
            "that way.")
    if abs(celda[0, 2]) > 1e-6 or abs(celda[1, 2]) > 1e-6:
        problemas.append(
            "the in-plane vectors have a z component. The cell has to "
            "be tetragonal\nor orthorhombic with z separated.")
    z = atoms.get_positions()[:, 2]
    largo = float(np.linalg.norm(celda[2]))
    if largo <= 0:
        problemas.append("the cell has no extent along z.")
    elif (z.max() - z.min()) > 0.98 * largo:
        problemas.append(
            "the atoms fill the whole cell along z. For a periodic "
            "electrode that is\nfine; for a scattering region there must be "
            "room left for the electrodes.")
    return problemas


def build_cond_input(prefixl: str, ikind: int = 0, prefixs: str = None,
                     prefixr: str = None, outdir: str = "./out",
                     energia0: float = 3.0, denergia: float = -0.1,
                     nenergia: int = 61, bdl: float = None,
                     bds: tuple = None, bdr: tuple = None,
                     band_file: str = "bands", tran_file: str = "trans.dat",
                     kpuntos=((0.0, 0.0, 1.0),), ewind: float = 1.0,
                     epsproj: float = 1e-3, nz1: int = 3,
                     ecut2d: float = None) -> str:
    lineas = [" &inputcond", f"    outdir='{outdir}'",
              f"    prefixl='{prefixl}'"]
    if prefixs:
        lineas.append(f"    prefixs='{prefixs}'")
    if prefixr:
        lineas.append(f"    prefixr='{prefixr}'")
    lineas += [f"    band_file='{band_file}'"]
    if ikind > 0:
        lineas.append(f"    tran_file='{tran_file}'")
    lineas += [f"    ikind={ikind}",
               f"    energy0={energia0}d0",
               f"    denergy={denergia}d0",
               f"    ewind={ewind}d0",
               f"    epsproj={epsproj:.1e}".replace("e-0", "d-0"),
               f"    nz1={nz1}"]
    if bdl is not None:
        lineas.append(f"    bdl={bdl:.6f}")
    if bds:
        lineas.append(f"    bds={bds[0]:.6f}")
    if bdr:
        lineas.append(f"    bdr={bdr[0]:.6f}")
    if ecut2d:
        lineas.append(f"    ecut2d={ecut2d}")
    lineas.append(" /")
    lineas.append(f"    {len(kpuntos)}")
    for k in kpuntos:
        lineas.append(f"    {k[0]} {k[1]} {k[2]}")
    lineas.append(f"    {nenergia}")
    return "\n".join(lineas) + "\n"


def prepare(electrodo, outdir: str = "balistico", dispersor=None,
            ikind: int = None, emin: float = -3.0, emax: float = 3.0,
            npuntos: int = 61, pseudo_dir: str = None,
            ecutwfc: float = None, ecutrho: float = None,
            kspacing: float = None, kpuntos=None,
            nz1: int = 3) -> tuple:
    """Escribe los scf de las regiones y el input de pwcond.x."""
    from qekit.modules import inputgen

    problemas = comprobar_geometria(electrodo)
    if dispersor is not None:
        problemas += comprobar_geometria(dispersor)
        c1 = np.array(electrodo.get_cell())[:2, :2]
        c2 = np.array(dispersor.get_cell())[:2, :2]
        if not np.allclose(c1, c2, atol=1e-4):
            problemas.append(
                "the electrode and the scattering region do NOT have the same "
                "cell in the xy plane.\npwcond.x joins the two regions "
                "there: if they do not match, there is no junction.")
    if problemas:
        raise ErrorDeUso("the geometry is not usable by pwcond.x:\n\n" +
                         "\n\n".join("  " + p for p in problemas))

    if ikind is None:
        ikind = 1 if dispersor is not None else 0
    if ikind == 2:
        raise ErrorDeUso(
            "ikind=2 (DIFFERENT left and right electrodes) is not "
            "implemented: Olla-DFT\nonly prepares the case of identical electrodes "
            "(ikind=1). For two different electrodes\nyou have to "
            "write the third scf and 'prefixr' and 'bdr' in cond.in by hand.")
    if ikind == 1 and dispersor is None:
        raise ErrorDeUso(
            "ikind=1 requires a scattering region: pass it with --scatterer, "
            "or use ikind=0\nto see only the complex bands of the electrode.")

    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    regiones = [("electrodo", electrodo)]
    if dispersor is not None:
        regiones.append(("dispersor", dispersor))

    prefijos = {}
    for nombre, at in regiones:
        common = sweep.prepare_common(at, pseudo_dir, ecutwfc, ecutrho,
                                      insulator=False, prefix=nombre[:6])
        grid = sweep.default_grid(at, kspacing)
        txt = inputgen.build_pw_input(
            atoms=at, pseudos=common["pseudos"], calculation="scf",
            prefix=common["prefix"], pseudo_dir=common["pseudo_dir"],
            ecutwfc=common["ecutwfc"], ecutrho=common["ecutrho"],
            kcard=f"K_POINTS automatic\n  {grid[0]} {grid[1]} {grid[2]} "
                  "0 0 0\n",
            insulator=False, degauss=common["degauss"],
            smearing=common["smearing"])
        sweep.write_input(out / f"scf_{nombre}.in", txt)
        prefijos[nombre] = common["prefix"]

    # el barrido de energia se hace de arriba abajo, como en los ejemplos
    paso = -(emax - emin) / max(npuntos - 1, 1)
    cond = build_cond_input(
        prefixl=prefijos["electrodo"], ikind=ikind,
        prefixs=prefijos.get("dispersor"),
        energia0=emax, denergia=paso, nenergia=npuntos,
        # fronteras derechas de cada región = longitud de SU celda en z
        bdl=None if dispersor is None else longitud_z(electrodo),
        bds=None if dispersor is None else (longitud_z(dispersor),),
        kpuntos=kpuntos or ((0.0, 0.0, 1.0),), nz1=nz1)
    sweep.write_input(out / "cond.in", cond)

    rep = ["--- Ballistic transport (pwcond.x) ---",
           f"Electrode: {electrodo.get_chemical_formula()} "
           f"({len(electrodo)} atoms)"]
    if dispersor is not None:
        rep.append(f"Scattering region: "
                   f"{dispersor.get_chemical_formula()} "
                   f"({len(dispersor)} atoms)")
    rep += [f"Mode: ikind = {ikind}"
            + ("  (complex bands only: the number of channels, which is the "
               "upper bound)" if ikind == 0 else
               "  (conductance with scattering region)"),
            f"Window: {emin} to {emax} eV relative to E_F, "
            f"{npuntos} points",
            "",
            f"Files in '{out.resolve()}':"]
    for nombre, _ in regiones:
        rep.append(f"  scf_{nombre}.in")
    rep += ["  cond.in", "",
            "Order:  " + "  &&  ".join(
                f"pw.x -in scf_{n}.in" for n, _ in regiones)
            + "  &&  pwcond.x -in cond.in",
            ""]
    if ikind == 0:
        rep += ["With ikind=0 the conductance is NOT produced: what comes out is the number of "
                "open channels at\neach energy, which is how much it COULD "
                "transmit at most. It is cheap and it is\nthe first thing to "
                "look at: if there are two channels at E_F, the conductance\ncannot "
                "exceed 2 G0.", ""]
    rep += [f"The conductance comes out in units of G0 = 2e²/h = "
            f"{G0:.4e} S,\nwhich corresponds to a resistance of "
            f"{R0 / 1000:.3f} kΩ per perfect channel.",
            "",
            "This is BALLISTIC transport: it applies to a nanocontact or a "
            "molecule between\nelectrodes, not to a macroscopic crystal. "
            "For that, use 'olla-dft transport'."]
    return {"prefijos": prefijos}, "\n".join(rep)


# ----------------------------------------------------------------------
# Lectura
# ----------------------------------------------------------------------


def collect(path) -> CondRun:
    p = Path(path)
    run = CondRun()

    trans = sorted(p.glob("trans*.dat")) + sorted(p.glob("*.tran"))
    if trans:
        d = np.loadtxt(trans[0], comments="#")
        if d.ndim == 1:
            d = d.reshape(1, -1)
        run.energias, run.transmision = d[:, 0], d[:, 1]
        run.ikind = 1

    salida = sorted(p.glob("cond*.out")) + sorted(p.glob("*.cond.out"))
    if salida:
        texto = salida[0].read_text(errors="ignore")
        m = re.search(r"ikind\s*=\s*(\d+)", texto)
        if m:
            run.ikind = int(m.group(1))
        if run.transmision is None:
            filas = re.findall(r"T_tot\s+(-?[\d.]+)\s+([\dEe.+-]+)", texto)
            if filas:
                d = np.array([[float(a), float(b)] for a, b in filas])
                run.energias, run.transmision = d[:, 0], d[:, 1]
        # canales abiertos: "Nchannels of the left tip = N", una vez por
        # energia y punto k. Si hay varios k por energia se toma el maximo.
        pares = re.findall(
            r"---\s+E-Ef\s*=\s*(-?[\d.]+).*?Nchannels of the left tip\s*=\s*"
            r"(\d+)", texto, re.S)
        if pares:
            d = np.array([[float(a), float(b)] for a, b in pares])
            energias = np.unique(d[:, 0])
            canales = np.array([d[d[:, 0] == e, 1].max() for e in energias])
            orden = np.argsort(energias)
            if run.energias is None:
                run.energias = energias[orden]
            run.canales = (canales[orden]
                           if len(energias) == len(run.energias) else None)

    if run.ikind is None:
        # pwcond.x no repite ikind en su salida; esta en el input, que
        # normalmente esta al lado. Sin el, el reporte no sabria si lo que
        # tiene delante es una conductancia o solo canales abiertos.
        entrada = p / "cond.in"
        if entrada.exists():
            m = re.search(r"ikind\s*=\s*(\d+)",
                          entrada.read_text(errors="ignore"))
            if m:
                run.ikind = int(m.group(1))

    if run.energias is None:
        raise ErrorDeUso(
            f"no pwcond.x result could be read in {p}.\n"
            "Look for trans*.dat (conductance) or the cond.out output (complex "
            "bands).")

    if run.transmision is not None:
        i = int(np.argmin(np.abs(run.energias)))
        run.G_fermi = float(run.transmision[i])
    _avisar(run)
    return run


def _avisar(run: CondRun) -> None:
    if run.transmision is not None and run.canales is not None and \
            len(run.canales) == len(run.transmision):
        exceso = run.transmision - run.canales
        if np.any(exceso > 0.01):
            run.avisos.append(
                "The transmission exceeds the number of open channels at "
                "some energy.\nThat is impossible: T <= N by construction. "
                "Check that the bdl/bds limits\nknow where each region "
                "ends and that the in-plane cells match.")
    if run.transmision is not None and np.any(run.transmision < -1e-6):
        run.avisos.append(
            "There are NEGATIVE transmissions. It is a sign that the calculation did not "
            "converge or that\nthe geometry of the regions is cut "
            "wrongly.")
    if run.ikind == 0:
        run.avisos.append(
            "Mode ikind=0: this is NOT the conductance. It is the number of "
            "open channels,\nwhich bounds the conductance from above. For "
            "the real value the scattering region\nand ikind=1 or 2 are needed.")


def report(run: CondRun) -> str:
    lines = ["--- Ballistic transport ---",
             f"Mode: ikind = {run.ikind}"]
    if run.energias is not None:
        lines.append(f"Window: {run.energias.min():.2f} to "
                     f"{run.energias.max():.2f} eV relative to E_F "
                     f"({len(run.energias)} points)")

    if run.G_fermi is not None:
        lines += ["",
                  "Conductance at the Fermi level:",
                  f"  T(E_F)  = {run.G_fermi:.4f}",
                  f"  G       = {run.G_fermi:.4f} G0 = "
                  f"{run.G_siemens:.4e} S",
                  f"  R       = {run.R_ohm / 1000:.3f} kΩ"]
        cerca = round(run.G_fermi)
        if cerca >= 1 and abs(run.G_fermi - cerca) < 0.08:
            lines.append(
                f"  T is very close to the integer {cerca}: that is {cerca} channel(s) "
                "transmitting almost\n  perfectly. This is conductance "
                "quantization, and seeing it is the\n  best "
                "sign that the calculation is well posed.")
        elif run.G_fermi < 0.1:
            lines.append(
                "  T almost zero: the contact is closed at that energy. With "
                "a molecule in\n  between this is normal if its gap falls on "
                "E_F.")

    if run.canales is not None:
        lines += ["", "Open channels (upper bound of T):",
                  f"  {'E-Ef (eV)':>11s} {'channels':>8s}"
                  + ("  {:>8s}".format("T") if run.transmision is not None
                     else "")]
        n = len(run.energias)
        idx = range(n) if n <= 12 else \
            [int(round(x)) for x in np.linspace(0, n - 1, 12)]
        for i in idx:
            fila = f"  {run.energias[i]:11.3f} {run.canales[i]:8.0f}"
            if run.transmision is not None and i < len(run.transmision):
                fila += f"  {run.transmision[i]:8.4f}"
            lines.append(fila)

    for a in run.avisos:
        lines += ["", a]
    lines += ["",
              f"G0 = 2e²/h = {G0:.4e} S; one perfect channel is "
              f"{R0 / 1000:.3f} kΩ.",
              "BALLISTIC transport: the electron crosses without scattering. "
              "For a macroscopic\ncrystal, which is diffusive, use "
              "'olla-dft transport'."]
    return "\n".join(lines)


def export(run: CondRun, outdir: str = ".") -> list:
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    f = out / "BALISTICO.dat"
    cols, nombres = [run.energias], ["E-EF(eV)"]
    if run.transmision is not None:
        cols.append(run.transmision); nombres.append("T")
    if run.canales is not None and len(run.canales) == len(run.energias):
        cols.append(run.canales); nombres.append("channels")
    np.savetxt(f, np.column_stack(cols), fmt="%14.6f",
               header=provenance.header_plain(
                   "ballistic transport",
                   {"ikind": run.ikind,
                    "G_en_G0": None if run.G_fermi is None
                    else round(run.G_fermi, 5)},
                   titulo="Landauer conductance") + "\n" +
               "  ".join(f"{n:>14s}" for n in nombres), comments="# ")
    txt = out / "BALISTICO.txt"
    txt.write_text(report(run) + "\n", encoding="utf-8")
    return [str(f), str(txt)]


def plot(run: CondRun, outfile: str = "balistico", formats="pdf,png",
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
    fig, ax = qstyle.new_figure(width, journal, 0.75)
    colores = qstyle.palette(4, mono=mono)

    if run.canales is not None and len(run.canales) == len(run.energias):
        ax.step(run.energias, run.canales, where="mid", lw=1.0,
                color=qstyle.INK_FAINT, dashes=[3, 2],
                label="open channels")
    if run.transmision is not None:
        ax.plot(run.energias, run.transmision, lw=1.5, color=colores[0],
                label="T(E)")
    ax.axvline(0.0, color=qstyle.INK_FAINT, lw=st["axis_line"],
               dashes=[3.5, 2.0])
    ax.set_xlabel(r"$E - E_\mathrm{F}$ (eV)")
    ax.set_ylabel(r"$T$  ($G/G_0$)")
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False, fontsize="small")
    return qstyle.save(fig, outfile, formats, dpi=dpi)
