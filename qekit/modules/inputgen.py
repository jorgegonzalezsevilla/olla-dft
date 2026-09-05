# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Generador de archivos de entrada para Quantum ESPRESSO (pw.x y post-proceso).

Presets disponibles:
  scf       -> scf.in
  relax     -> relax.in            (posiciones atómicas)
  vc-relax  -> vc-relax.in         (posiciones + celda)
  nscf      -> nscf.in             (malla densa, tetraedros)
  bands     -> scf.in + bands.in + bands_pp.in + KPATH.txt   (celda primitiva)
  dos       -> scf.in + nscf.in + dos.in + projwfc.in
  all       -> flujo completo: scf, nscf, bands y post-proceso
  md        -> md.in               (dinámica molecular Born-Oppenheimer)
Todos generan además run.sh y run.py con el orden de ejecución.
"""

import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.data import atomic_masses, atomic_numbers

from qekit import config as qcfg
from qekit.core import kpoints, plataforma, pseudo, structure
from qekit.core.errors import ErrorDeUso

PRESETS = ("scf", "relax", "vc-relax", "nscf", "bands", "dos", "all", "md")

# 1 unidad de tiempo atómica de Rydberg = 4.8378e-2 fs. pw.x pide dt en esas
# unidades y en cambio IMPRIME el paso en fs, que es lo que lee 'olla-dft md';
# la conversión vive aquí para que las dos mitades no se desincronicen.
_FS_POR_UA = 4.8378e-2

# Nombres de vdw_corr que entiende pw.x. 'DFT-D' es el alias viejo de
# grimme-d2 y se acepta porque aparece así en inputs publicados.
VDW = ("grimme-d2", "grimme-d3", "DFT-D", "ts-vdw", "xdm", "mbd")

# Funcionales híbridos que entiende pw.x, con la fracción de intercambio
# exacto y el parámetro de apantallamiento de cada uno. QE trae los valores
# por defecto, pero escribirlos explícitamente deja el input auto-explicativo
# y permite cambiarlos sin tocar el código.
HIBRIDOS = {
    "hse":   {"input_dft": "hse",   "exx_fraction": 0.25,
              "screening_parameter": 0.106,
              "nombre": "HSE06 (screened, α=0.25, ω=0.106 bohr⁻¹)"},
    "pbe0":  {"input_dft": "pbe0",  "exx_fraction": 0.25,
              "nombre": "PBE0 (α=0.25, unscreened)"},
    "b3lyp": {"input_dft": "b3lyp", "exx_fraction": 0.20,
              "nombre": "B3LYP (α=0.20; intended for molecules)"},
    "gaupbe": {"input_dft": "gaupbe", "exx_fraction": 0.24,
               "nombre": "Gau-PBE (screened with a Gaussian)"},
}


def espesor_celda(atoms, eje: int) -> float:
    """Altura de la celda en Å a lo largo de la normal al plano de los otros
    dos vectores. Para una celda ortogonal es simplemente |c|, pero en una
    celda inclinada es lo único que mide de verdad cuánto vacío hay."""
    import numpy as np
    cell = np.asarray(atoms.cell.array, dtype=float)
    otros = [i for i in range(3) if i != eje]
    area = np.linalg.norm(np.cross(cell[otros[0]], cell[otros[1]]))
    vol = abs(np.linalg.det(cell))
    return float(vol / area) if area > 0 else 0.0


def hueco_vacio(atoms, eje: int):
    """Hueco más ancho entre átomos a lo largo de `eje` (0, 1, 2).

    Devuelve (centro_fraccionario, hueco_fraccionario, hueco_en_angstrom).
    Se mide en Å y no en fraccionarias a propósito: en una celda primitiva
    inclinada dos átomos vecinos por la diagonal dejan un hueco fraccionario
    enorme a lo largo de un eje sin que exista vacío ninguno.
    """
    d = espesor_celda(atoms, eje)
    frac = sorted(float(f) % 1.0 for f in atoms.get_scaled_positions()[:, eje])
    if not frac:
        return 0.5, 1.0, d
    if len(frac) == 1:
        return (frac[0] + 0.5) % 1.0, 1.0, d
    huecos = [(frac[i + 1] - frac[i], frac[i]) for i in range(len(frac) - 1)]
    huecos.append((1.0 - frac[-1] + frac[0], frac[-1]))
    hueco, inicio = max(huecos)
    return (inicio + hueco / 2.0) % 1.0, hueco, hueco * d


def _divisor(n: int) -> int:
    """El mayor divisor propio de n, o 1. Para sugerir mallas de EXX válidas."""
    for d in range(n // 2, 0, -1):
        if n % d == 0:
            return d
    return 1


def _region_vacio(atoms, edir: int):
    """Dónde poner el diente de sierra de la corrección dipolar.

    Devuelve (emaxpos, eopreg) en coordenadas fraccionarias del eje `edir`.
    La sierra tiene que subir y bajar DENTRO del vacío: se pone el máximo en
    el centro del hueco más ancho y se le da a la bajada un tercio del hueco.
    """
    eje = int(edir) - 1
    centro, hueco, hueco_A = hueco_vacio(atoms, eje)
    if hueco_A < 5.0:
        raise ErrorDeUso(
            "the dipole correction needs vacuum along direction "
            f"{'abc'[eje]}, and the widest gap between atoms is "
            f"{hueco_A:.1f} Å. If this is a slab, add vacuum (for example "
            "'olla-dft surface --vacuum 20'); if it is bulk material, there is no "
            "dipole to correct and --dipole is unnecessary. With the sawtooth crossing "
            "the atoms the result is worse than without correction.")
    return centro, min(0.2, max(0.02, hueco / 3.0))


# ----------------------------------------------------------------------
# Opciones del generador
# ----------------------------------------------------------------------
@dataclass
class GenOptions:
    preset: str = "scf"
    outdir: str = "."            # carpeta donde se escriben los archivos
    kspacing: float = 0.20       # malla scf (Å^-1, incluye 2π)
    kgrid: tuple = None          # malla scf explícita (n1, n2, n3); anula kspacing
    kspacing_nscf: float = 0.12  # malla nscf/DOS
    band_points: int = 20        # puntos por segmento del k-path
    ecutwfc: float = None        # None => automático (UPF o defaults)
    ecutrho: float = None
    insulator: bool = False      # True => occupations='fixed' en scf
    use_primitive: bool = None   # None => automático (bands/all sí)
    pseudo_dir: str = None       # None => tomado de la configuración
    prefix: str = None           # None => fórmula química reducida
    nspin: int = 1               # 2 => cálculo con polarización de espín
    magnetization: dict = None   # símbolo -> magnetización inicial (fracción)
    vdw: str = None              # 'grimme-d2', 'grimme-d3', 'xdm', 'ts-vdw'
    soc: bool = False            # acoplamiento espín-órbita (pseudos rel.)
    hubbard: dict = None         # símbolo -> U en eV
    hubbard_style: str = "legacy"
    tot_charge: float = None     # carga total de la celda (electrones de más)
    dipole: object = False       # False, True (eje c) o 1/2/3
    nosym: bool = False          # desactivar simetría
    md: dict = None              # dt_fs, nstep, thermostat, temperature
    hibrido: str = None          # 'hse', 'pbe0', 'b3lyp', 'gaupbe'
    exx_grid: tuple = None       # malla q del intercambio exacto
    exx_fraction: float = None
    notes: list = field(default_factory=list)


def parse_magnetization(text: str, symbols: list) -> dict:
    """Interpreta el argumento --mag.

    Acepta un solo número ('0.5', aplicado a todos los elementos) o pares
    por elemento ('Fe=0.7,O=0'). Devuelve símbolo -> valor.
    """
    if not text:
        return {}
    text = text.strip()
    unique = list(dict.fromkeys(symbols))
    if "=" not in text:
        try:
            value = float(text)
        except ValueError:
            raise ErrorDeUso(
                f"could not parse the magnetization '{text}'. "
                "Use a number (0.5) or per-element pairs (Fe=0.7,O=0)."
            )
        return {s: value for s in unique}
    mag = {}
    for chunk in text.replace(";", ",").split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        if "=" not in chunk:
            raise ErrorDeUso(f"invalid format in '{chunk}'; expected El=value.")
        el, val = chunk.split("=", 1)
        el = el.strip()
        if el not in unique:
            raise ErrorDeUso(
                f"element '{el}' is not in the structure "
                f"({', '.join(unique)})."
            )
        try:
            mag[el] = float(val)
        except ValueError:
            raise ErrorDeUso(
                f"could not parse the magnetization '{val}' for {el}. "
                "Use a number (for example Fe=0.7,O=0)."
            ) from None
    return mag


# ----------------------------------------------------------------------
# Utilidades de formato Fortran
# ----------------------------------------------------------------------
def _fval(value) -> str:
    if isinstance(value, bool):
        return ".true." if value else ".false."
    if isinstance(value, float):
        return f"{value:g}"
    if isinstance(value, int):
        return str(value)
    return f"'{value}'"


def _namelist(name: str, params: dict) -> str:
    lines = [f"&{name}"]
    for key, value in params.items():
        if value is None:
            continue
        lines.append(f"  {key:16s} = {_fval(value)}")
    lines.append("/")
    return "\n".join(lines) + "\n"


# ----------------------------------------------------------------------
# Tarjetas de estructura
# ----------------------------------------------------------------------
def fixed_atoms(atoms: Atoms) -> set:
    """Índices de los átomos que no deben moverse en un relax.

    Se aceptan las dos formas en que puede llegar la marca: el array
    `qekit_fijo` que deja 'surface --fix' y la restricción `FixAtoms` de
    ASE (que sobrevive a un POSCAR o a un input de pw.x, cosa que el CIF no
    hace). Devuelve un conjunto de índices, vacío si no hay nada fijo.
    """
    fijos = set()
    marca = atoms.arrays.get("qekit_fijo")
    if marca is not None:
        fijos.update(int(i) for i in np.flatnonzero(np.asarray(marca) != 0))
    for c in atoms.constraints or []:
        idx = getattr(c, "index", None)
        if idx is None or c.__class__.__name__ != "FixAtoms":
            continue
        idx = np.asarray(idx)
        if idx.dtype == bool:
            idx = np.flatnonzero(idx)
        fijos.update(int(i) for i in idx)
    return fijos


def _structure_cards(atoms: Atoms, pseudos: dict,
                     extra_species: list = None) -> str:
    lines = ["ATOMIC_SPECIES"]
    for symbol in dict.fromkeys(atoms.get_chemical_symbols()):
        mass = atomic_masses[atomic_numbers[symbol]]
        lines.append(f"  {symbol:3s} {mass:10.4f}  {pseudos[symbol]['filename']}")
    # Especies declaradas SIN átomos que las usen. Suena raro, pero es
    # exactamente lo que pide initial_state.x: la contraparte con hueco de
    # core existe como TIPO para que excite(nt) pueda apuntarle, y ningún
    # átomo la ocupa.
    for etiqueta, elemento, upf in (extra_species or []):
        mass = atomic_masses[atomic_numbers[elemento]]
        lines.append(f"  {etiqueta:3s} {mass:10.4f}  {upf}")
    lines.append("")
    lines.append("ATOMIC_POSITIONS crystal")
    fijos = fixed_atoms(atoms)
    for i, (symbol, spos) in enumerate(zip(
        atoms.get_chemical_symbols(), atoms.get_scaled_positions()
    )):
        # tercera columna: 0 0 0 = el átomo no se mueve en relax/md
        lines.append(
            f"  {symbol:3s} {spos[0]:14.10f} {spos[1]:14.10f} {spos[2]:14.10f}"
            + ("   0 0 0" if i in fijos else "")
        )
    lines.append("")
    lines.append("CELL_PARAMETERS angstrom")
    for vec in atoms.cell.array:
        lines.append(f"  {vec[0]:14.10f} {vec[1]:14.10f} {vec[2]:14.10f}")
    return "\n".join(lines) + "\n"


def _kgrid_card(grid: tuple) -> str:
    if grid == (1, 1, 1):
        return "K_POINTS gamma\n"
    n1, n2, n3 = grid
    return f"K_POINTS automatic\n  {n1} {n2} {n3} 0 0 0\n"


# ----------------------------------------------------------------------
# Estimación de nbnd para nscf/bands
# ----------------------------------------------------------------------
def _estimate_nbnd(atoms: Atoms, pseudos: dict):
    total = 0.0
    for symbol in atoms.get_chemical_symbols():
        z = pseudos[symbol].get("z_valence")
        if z is None:
            return None  # sin UPF no se puede estimar
        total += z
    occupied = total / 2.0
    return int(math.ceil(occupied * 1.25 + 4))


# ----------------------------------------------------------------------
# Construcción de un input de pw.x
# ----------------------------------------------------------------------
def build_pw_input(
    atoms: Atoms,
    pseudos: dict,
    calculation: str,
    prefix: str,
    pseudo_dir: str,
    ecutwfc: float,
    ecutrho: float,
    kcard: str,
    insulator: bool = False,
    degauss: float = 0.01,
    smearing: str = "cold",
    nbnd: int = None,
    tetrahedra: bool = False,
    nspin: int = 1,
    magnetization: dict = None,
    tot_magnetization: float = None,
    conv_thr: float = 1.0e-8,
    forc_conv_thr: float = 1.0e-4,
    etot_conv_thr: float = 1.0e-5,
    vdw: str = None,
    nosym: bool = False,
    disk_io: str = None,
    soc: bool = False,
    hubbard: dict = None,
    hubbard_style: str = 'legacy',
    tot_charge: float = None,
    dipole_correction=False,
    md: dict = None,
    cell_dofree: str = None,
    hibrido: str = None,
    exx_grid: tuple = None,
    exx_fraction: float = None,
    extra_species: list = None,
) -> str:
    nat = len(atoms)
    species = list(dict.fromkeys(atoms.get_chemical_symbols()))
    ntyp = len(species) + len(extra_species or [])

    control = {
        "calculation": calculation,
        "prefix": prefix,
        "outdir": "./out",
        "pseudo_dir": pseudo_dir,
        "tprnfor": True,
        "tstress": True,
    }
    if disk_io:
        # 'nowf' evita escribir las funciones de onda. Para un nscf de malla
        # densa del que solo se quieren los AUTOVALORES (transporte,
        # superficie de Fermi) eso ahorra cientos de MB y bastante tiempo.
        control["disk_io"] = disk_io
    if calculation in ("relax", "vc-relax"):
        control["etot_conv_thr"] = etot_conv_thr
        control["forc_conv_thr"] = forc_conv_thr
    if calculation in ("md", "vc-md"):
        md = dict(md or {})
        # dt de pw.x va en unidades de Rydberg atómicas (1 u.a. = 4.8378e-2 fs)
        control["dt"] = round(float(md.get("dt_fs", 1.0)) / _FS_POR_UA, 4)
        control["nstep"] = int(md.get("nstep", 1000))

    system = {
        "ibrav": 0,
        "nat": nat,
        "ntyp": ntyp,
        "ecutwfc": ecutwfc,
        "ecutrho": ecutrho,
    }
    if soc:
        # SOC en QE = cálculo no colineal + lspinorb. Exige pseudos
        # totalmente relativistas; quien llama debe haberlo verificado.
        system["noncolin"] = True
        system["lspinorb"] = True
        system.pop("nspin", None)
    if hibrido:
        h = HIBRIDOS.get(str(hibrido).lower())
        if h is None:
            raise ErrorDeUso(
                f"unknown hybrid functional '{hibrido}'. "
                f"Options: {', '.join(sorted(HIBRIDOS))}.")
        system["input_dft"] = h["input_dft"]
        system["exx_fraction"] = (float(exx_fraction) if exx_fraction is not None
                                  else h["exx_fraction"])
        if "screening_parameter" in h:
            system["screening_parameter"] = h["screening_parameter"]
        # La malla de q para el intercambio exacto. Es LO que decide el coste:
        # cada punto q multiplica el trabajo, y por eso lo normal es usarla
        # más gruesa que la de k. Tiene que DIVIDIR la malla de k, o pw.x se
        # queja de "nqx must be a divisor of nk".
        q = tuple(int(x) for x in (exx_grid or (1, 1, 1)))
        system["nqx1"], system["nqx2"], system["nqx3"] = q
    if tot_charge is not None:
        system["tot_charge"] = float(tot_charge)
    if dipole_correction:
        # Losa polar: el dipolo neto de la losa se ve a sí mismo a través del
        # vacío por las condiciones periódicas y desplaza las dos funciones
        # trabajo hacia un promedio falso. QE lo cancela con un potencial en
        # diente de sierra, pero SOLO si la sierra cae dentro del vacío: si
        # cruza los átomos, la corrección es peor que no ponerla.
        edir = 3 if dipole_correction is True else int(dipole_correction)
        if edir not in (1, 2, 3):
            raise ErrorDeUso(
                f"the dipole-correction direction is 1, 2 or 3 (a, b or c); "
                f"got {edir}.")
        emaxpos, eopreg = _region_vacio(atoms, edir)
        system["dipfield"] = True
        system["edir"] = edir
        system["emaxpos"] = round(emaxpos, 4)
        system["eopreg"] = round(eopreg, 4)
        system["eamp"] = 0.0          # solo corregir, sin campo aplicado
        control["tefield"] = True
    if hubbard and hubbard_style == "legacy":
        # sintaxis de QE <= 7.0: lda_plus_u + Hubbard_U(i) en &SYSTEM,
        # indexado por el ORDEN de ATOMIC_SPECIES
        system["lda_plus_u"] = True
        for i, sym in enumerate(species, start=1):
            if sym in hubbard:
                system[f"Hubbard_U({i})"] = float(hubbard[sym])
    if nosym:
        # epsilon.x y el transporte necesitan la malla COMPLETA de k: con la
        # reducción por simetría los puntos no forman una rejilla y ni la
        # suma sobre la zona ni las derivadas de E(k) salen bien.
        system["nosym"] = True
        system["noinv"] = True
    if nbnd:
        system["nbnd"] = nbnd
    if nspin == 2:
        system["nspin"] = 2
        if tot_magnetization is not None:
            # Con occupations='fixed' y espín, QE exige saber cuántos
            # electrones van a cada canal: sin esto aborta con "the system is
            # metallic, specify occupations" incluso en un aislante, porque
            # con un número impar de electrones no puede repartirlos solo.
            system["tot_magnetization"] = float(tot_magnetization)
        mag = magnetization or {}
        for i, sym in enumerate(species, start=1):
            # QE indexa starting_magnetization por el orden de ATOMIC_SPECIES
            system[f"starting_magnetization({i})"] = float(mag.get(sym, 0.0))
    if vdw:
        # Nombres que entiende QE: 'grimme-d2' (alias DFT-D, desde QE 5),
        # 'grimme-d3' (desde QE 7.1). Se pasa tal cual.
        system["vdw_corr"] = vdw
    if tetrahedra:
        system["occupations"] = "tetrahedra_opt"
    elif insulator:
        system["occupations"] = "fixed"
    else:
        system["occupations"] = "smearing"
        system["smearing"] = smearing
        system["degauss"] = degauss

    # mixing_beta: 0.7 (valor por omisión de QE) converge en menos iteraciones
    # en aislantes con ocupaciones fijas; 0.4 sigue siendo la opción prudente
    # con smearing, donde la oscilación de carga es el riesgo real.
    electrons = {
        "conv_thr": conv_thr,
        "mixing_beta": 0.7 if (insulator and not tetrahedra) else 0.4,
        "electron_maxstep": 200,
    }

    text = _namelist("CONTROL", control) + _namelist("SYSTEM", system)
    text += _namelist("ELECTRONS", electrons)
    if calculation in ("relax", "vc-relax"):
        text += _namelist("IONS", {"ion_dynamics": "bfgs"})
    if calculation in ("md", "vc-md"):
        ions = {"ion_dynamics": "verlet",
                "pot_extrapolation": "second_order",
                "wfc_extrapolation": "second_order"}
        term = str(md.get("thermostat", "none")).lower()
        if term != "none":
            ions["ion_temperature"] = term
            ions["tempw"] = float(md.get("temperature", 300.0))
            if term == "rescaling":
                ions["tolp"] = float(md.get("tolp", 50.0))
            if term in ("berendsen", "andersen", "svr"):
                ions["nraise"] = int(md.get("nraise", 100))
        text += _namelist("IONS", ions)
    if calculation in ("vc-relax", "vc-md"):
        celda = {"cell_dynamics": "bfgs" if calculation == "vc-relax" else "pr",
                 "press_conv_thr": 0.05}
        if cell_dofree:
            # Relajar solo parte de la celda. En un barrido de deformación es
            # lo que permite fijar la deformación impuesta y dejar que el eje
            # perpendicular responda (relajación de Poisson): sin esto, una
            # monocapa estirada en el plano conserva un espesor que ya no le
            # corresponde y la energía sale alta por un motivo falso.
            celda["cell_dofree"] = cell_dofree
        text += _namelist("CELL", celda)
    text += "\n" + _structure_cards(atoms, pseudos, extra_species) + "\n"
    if hubbard and hubbard_style == "card":
        # sintaxis de QE >= 7.1: tarjeta HUBBARD con el orbital explícito
        text += "HUBBARD (ortho-atomic)\n"
        for sym, u in hubbard.items():
            text += f"  U {sym}-{_orbital_hubbard(sym)} {float(u):g}\n"
        text += "\n"
    text += kcard
    return text


# orbital que se corrige con U, por elemento (el de la capa d o f abierta)
def _orbital_hubbard(symbol: str) -> str:
    from ase.data import atomic_numbers
    z = atomic_numbers.get(symbol, 0)
    if 57 <= z <= 71 or 89 <= z <= 103:
        return "4f" if z <= 71 else "5f"
    if 21 <= z <= 30:
        return "3d"
    if 39 <= z <= 48:
        return "4d"
    if 72 <= z <= 80:
        return "5d"
    return "2p"


# ----------------------------------------------------------------------
# Inputs de post-proceso (dos.x, projwfc.x, bands.x)
# ----------------------------------------------------------------------
def build_dos_input(prefix: str) -> str:
    return _namelist(
        "DOS",
        {"prefix": prefix, "outdir": "./out", "fildos": f"{prefix}.dos", "DeltaE": 0.02},
    )


def build_projwfc_input(prefix: str) -> str:
    return _namelist(
        "PROJWFC",
        {"prefix": prefix, "outdir": "./out", "filpdos": f"{prefix}.pdos", "DeltaE": 0.02},
    )


def build_bandsx_input(prefix: str) -> str:
    return _namelist(
        "BANDS",
        {"prefix": prefix, "outdir": "./out", "filband": f"{prefix}.bands.dat", "lsym": True},
    )


# ----------------------------------------------------------------------
# Script de ejecución
# ----------------------------------------------------------------------
_RUN_CMDS = {
    "scf.in": "pw.x",
    "relax.in": "pw.x",
    "vc-relax.in": "pw.x",
    "md.in": "pw.x",
    "nscf.in": "pw.x",
    "bands.in": "pw.x",
    "bands_pp.in": "bands.x",
    "dos.in": "dos.x",
    "projwfc.in": "projwfc.x",
}


def build_run_script(input_files: list, nproc: int) -> str:
    lines = [
        "#!/bin/bash",
        "# Generated by Olla-DFT — runs the calculations in order.",
        "# On Windows, or without bash:  python run.py",
        # Sin pipefail, `pw.x | tee` devuelve el código de tee y un fallo de
        # QE parece un éxito. Bash es deliberado: POSIX sh no tiene pipefail.
        "set -e -o pipefail",
        f"NP=${{NPROC:-{nproc}}}",
        'if [ "$NP" -gt 1 ] && command -v mpirun >/dev/null 2>&1; then',
        '  LANZ="mpirun -np $NP"',
        'elif [ "$NP" -gt 1 ] && command -v mpiexec >/dev/null 2>&1; then',
        '  LANZ="mpiexec -n $NP"',
        "else",
        '  LANZ=""',
        "fi",
        "",
    ]
    for fname in input_files:
        exe = _RUN_CMDS.get(fname)
        if exe is None:
            continue
        out = fname.rsplit(".", 1)[0] + ".out"
        lines.append(f'echo ">> {exe} < {fname}"')
        lines.append(f'$LANZ {exe} -in {fname} | tee {out}')
        lines.append("")
    return "\n".join(lines) + "\n"


def build_run_python_script(input_files: list, nproc: int) -> str:
    """Versión portable de ``run.sh`` para Windows o máquinas sin Bash."""
    pasos = []
    for fname in input_files:
        exe = _RUN_CMDS.get(fname)
        if exe is None:
            continue
        out = fname.rsplit(".", 1)[0] + ".out"
        pasos.append((exe, fname, out))
    return plataforma.build_sequential_python_script(pasos, nproc)


# ----------------------------------------------------------------------
# Orquestador
# ----------------------------------------------------------------------
def generate(atoms: Atoms, opts: GenOptions) -> str:
    """Genera todos los archivos del preset y devuelve un reporte legible."""
    # Validate before creating directories or writing any input.
    if opts.kgrid is not None:
        from numbers import Integral
        if (len(opts.kgrid) != 3 or any(isinstance(n, bool) or
                not isinstance(n, Integral) or n < 1 for n in opts.kgrid)):
            raise ErrorDeUso("--kgrid needs three positive integers")
    for name in ("ecutwfc", "ecutrho"):
        value = getattr(opts, name)
        if value is not None and (not np.isfinite(value) or value <= 0):
            raise ErrorDeUso(f"--{name} must be positive and finite")
    for name in ("kspacing", "kspacing_nscf"):
        value = getattr(opts, name)
        if value is not None and (not np.isfinite(value) or value < 0):
            raise ErrorDeUso(f"--{name} must be non-negative and finite")
    cfg = qcfg.load()
    preset = opts.preset
    if preset not in PRESETS:
        raise ErrorDeUso(f"unknown preset '{preset}'. Options: {', '.join(PRESETS)}")

    pseudo_dir = opts.pseudo_dir or cfg["pseudo_dir"]
    degauss = float(cfg["degauss"])
    smearing = cfg["smearing"]
    nproc = int(cfg["nproc"])
    report = []

    needs_bands = preset in ("bands", "all")
    use_primitive = opts.use_primitive
    if use_primitive is None:
        use_primitive = needs_bands

    # --- celda de trabajo y k-path ---
    kpath = None
    if needs_bands:
        kpath = kpoints.get_kpath(atoms)
        work_atoms = kpath.primitive  # el k-path está referido a esta celda
        if kpath.cell_changed:
            report.append(
                "WARNING: the standardized primitive cell (seekpath) was used in all\n"
                "inputs, because the band k-path is referred to it."
            )
    elif use_primitive:
        work_atoms = structure.primitive(atoms)
        if len(work_atoms) != len(atoms):
            report.append(
                f"Cell reduced to the primitive: {len(atoms)} -> {len(work_atoms)} atoms."
            )
    else:
        work_atoms = atoms

    # --- pseudopotenciales y cutoffs ---
    symbols = work_atoms.get_chemical_symbols()
    # Pasa por el mismo selector que el resto de Olla-DFT: respeta los
    # --pseudo del usuario y evita que gane el primero por orden
    # alfabetico, que es como se cuela un O de BLYP junto a un Ni de PBE.
    from qekit.modules import sweep as _sweep
    pseudos = pseudo.resolve(symbols, pseudo_dir,
                             forzados=_sweep.pseudo_overrides())
    if opts.soc:
        # lspinorb con pseudos escalares no falla: da un desdoblamiento
        # espín-órbita de CERO que parece un resultado. Se corta aquí, antes
        # de escribir nada (es la misma comprobación que usa `sweep`).
        _sweep.check_soc_pseudos({"pseudos": pseudos})
    missing = [s for s, p in pseudos.items() if not p["found"]]
    default_wfc = float(cfg["ecutwfc"])
    dual = float(cfg["dual"])
    auto_wfc, auto_rho = pseudo.recommend_cutoffs(pseudos, default_wfc, dual)
    ecutwfc = opts.ecutwfc if opts.ecutwfc is not None else auto_wfc
    ecutrho = opts.ecutrho if opts.ecutrho is not None else auto_rho

    prefix = opts.prefix or work_atoms.get_chemical_formula(
        mode="hill", empirical=True
    )

    report.append(f"Structure: {work_atoms.get_chemical_formula()} "
                  f"({len(work_atoms)} atoms)  |  prefix = '{prefix}'")
    report.append(f"Pseudopotentials in: {pseudo_dir}")
    for sym, p in pseudos.items():
        if p["found"]:
            extra = ""
            if p["ecutwfc"]:
                extra = f"  (suggested cutoff: {p['ecutwfc']:.0f} Ry)"
            report.append(f"  {sym:3s} -> {p['filename']}{extra}")
            if p["alternatives"]:
                report.append(f"        other options: {', '.join(p['alternatives'])}")
        else:
            report.append(f"  {sym:3s} -> NOT FOUND ('{p['filename']}' was written)")
    if missing:
        report.append(
            "ATTENTION: pseudopotentials are missing. Download them (e.g. from SSSP,\n"
            "https://www.materialscloud.org/sssp) and place them in pseudo_dir,\n"
            "or adjust the path with:  olla-dft config set pseudo_dir /path/to/pseudos"
        )
    report.append(f"Cutoffs: ecutwfc = {ecutwfc} Ry, ecutrho = {ecutrho} Ry"
                  + ("  (automatic)" if not opts.ecutwfc else ""))

    # --- mallas de k-points ---
    grid_scf = tuple(opts.kgrid) if opts.kgrid else kpoints.kgrid_from_spacing(work_atoms, opts.kspacing)
    grid_nscf = kpoints.kgrid_from_spacing(work_atoms, opts.kspacing_nscf)
    report.append(f"k-mesh (scf):  {grid_scf[0]} x {grid_scf[1]} x {grid_scf[2]}")
    if preset in ("nscf", "dos", "all"):
        report.append(
            f"k-mesh (nscf): {grid_nscf[0]} x {grid_nscf[1]} x {grid_nscf[2]}"
        )

    nbnd = None
    if preset in ("nscf", "bands", "dos", "all"):
        nbnd = _estimate_nbnd(work_atoms, pseudos)
        if nbnd and opts.nspin == 2:
            # con espín conviene un margen mayor de bandas vacías
            nbnd = int(nbnd * 1.2) + 2
        if nbnd:
            report.append(f"Estimated nbnd for nscf/bands: {nbnd}")

    if opts.nspin == 2:
        mag_txt = ", ".join(
            f"{s}={v:g}" for s, v in (opts.magnetization or {}).items()
        ) or "0 (no initial magnetization)"
        report.append(f"Spin polarization: enabled  |  initial magnetization: {mag_txt}")
        if not opts.magnetization:
            report.append(
                "  WARNING: without an initial magnetization the calculation usually converges\n"
                "  to the non-magnetic solution. Use --mag (for example --mag Fe=0.7)."
            )

    if opts.vdw:
        report.append(f"van der Waals correction: vdw_corr = '{opts.vdw}'")
    if opts.soc:
        report.append("Spin-orbit: non-collinear calculation with lspinorb "
                      "(requires fully relativistic pseudopotentials)")
    if opts.hubbard:
        u_txt = ", ".join(f"{k}={v:g} eV" for k, v in opts.hubbard.items())
        report.append(f"Hubbard U: {u_txt}  (syntax {opts.hubbard_style})")
        report.append("  A hand-set U is a choice, not a result. "
                      "To compute it:  olla-dft hubbard --cycle")
    if opts.hibrido:
        h = HIBRIDOS.get(str(opts.hibrido).lower())
        if h is None:
            raise ErrorDeUso(
                f"unknown hybrid functional '{opts.hibrido}'. "
                f"Options: {', '.join(sorted(HIBRIDOS))}.")
        q = tuple(int(x) for x in (opts.exx_grid or (1, 1, 1)))
        malos = [(qi, ki, eje) for qi, ki, eje in zip(q, grid_scf, "abc")
                 if ki % qi]
        if malos:
            detalle = "; ".join(f"nqx{eje} = {qi} does not divide {ki}"
                                for qi, ki, eje in malos)
            raise ErrorDeUso(
                f"the exact-exchange mesh has to DIVIDE the k-mesh, and "
                f"{detalle}. pw.x stops with 'nqx must be a divisor of "
                f"nk'. With the k-mesh {grid_scf[0]}x{grid_scf[1]}x{grid_scf[2]} "
                f"valid choices are, for example, "
                f"{'x'.join(str(_divisor(k)) for k in grid_scf)} or 1x1x1.")
        nq = q[0] * q[1] * q[2]
        report.append(f"Hybrid functional: {h['nombre']}")
        report.append(f"  Exact-exchange mesh: {q[0]}x{q[1]}x{q[2]} "
                      f"= {nq} q-point{'s' if nq > 1 else ''}")
        # Coste MEDIDO, no supuesto. Silicio de 2 átomos, malla k 4x4x4,
        # 20 Ry, el mismo binario y la misma máquina:
        #     LDA          0.69 s
        #     HSE nq=1     3.98 s   (x5.7)
        #     HSE nq=8    16.2  s   (x23.5)
        #     HSE nq=64  116.4  s   (x168)
        # es decir, factor ≈ 3 + 2.6·n_q. En celdas mayores el factor CRECE,
        # porque el intercambio exacto va con el número de parejas de bandas
        # ocupadas y el término semilocal no.
        factor = 3.0 + 2.6 * nq
        report.append(
            f"  COST: about {factor:.0f} times the same calculation with PBE.\n"
            f"  Measured on 2-atom silicon: the factor is "
            f"approximately 3 + 2.6·n_q, and\n  GROWS with the size of the "
            f"cell, because exact exchange scales with the\n  pairs of "
            f"occupied bands. Estimate the PBE one (--estimate) and multiply.")
        report.append(
            "  CONVERGENCE: the q-mesh is not a numerical detail, it changes the "
            "result.\n  In that same silicon the gap came out 2.68 eV with "
            "nqx 1x1x1, 1.83 with 2x2x2 and\n  1.41 with 4x4x4 (the "
            "experimental value is 1.17). A coarse q-mesh does NOT give an\n  "
            "approximate hybrid: it gives a gap that is too large. Increase it until the "
            "number\n  stops moving, or do not cite the hybrid."
            + ("\n  With 1x1x1 the result will come out clearly overestimated: "
               "it serves to check that\n  the calculation runs, not for a number."
               if nq == 1 else ""))
        if preset in ("bands", "all"):
            report.append(
                "  WARNING: bands with a hybrid. pw.x can NOT do a "
                "'bands' calculation with EXX\n  from the density: you have to "
                "interpolate (Wannier90) or run an scf with\n  the k-points of the "
                "path included. This input will give you the scf, not the "
                "dispersion.")
    if opts.tot_charge is not None:
        q = float(opts.tot_charge)
        report.append(
            f"Total charge: tot_charge = {q:g} "
            f"({'missing' if q > 0 else 'excess'} {abs(q):g} electrons)")
        report.append("  QE compensates with a uniform background: the energy of a "
                      "charged cell is NOT comparable with the neutral one without correcting "
                      "the finite size.")
    if opts.dipole:
        eje = 3 if opts.dipole is True else int(opts.dipole)
        emaxpos, eopreg = _region_vacio(work_atoms, eje)
        report.append(
            f"Dipole correction: edir = {eje} ({'abc'[eje - 1]}), "
            f"sawtooth at emaxpos = {emaxpos:.3f}, eopreg = {eopreg:.3f}")
    if preset == "md":
        md = dict(opts.md or {})
        dt_fs = float(md.get("dt_fs", 1.0))
        nstep = int(md.get("nstep", 1000))
        term = str(md.get("thermostat", "none")).lower()
        report.append(
            f"Molecular dynamics: {nstep} steps of {dt_fs:g} fs "
            f"= {nstep * dt_fs / 1000.0:.2f} ps"
            + (f", thermostat '{term}' at {float(md.get('temperature', 300.0)):g} K"
               if term != "none" else ", microcanonical (NVE, no thermostat)"))
        report.append("  Symmetry disabled (nosym): mandatory in MD.")
        if len(work_atoms) < 20:
            report.append(
                f"  WARNING: {len(work_atoms)} atoms are few for an MD run. "
                "With a small cell the periodic replicas move in "
                "phase and g(r) and the MSD come out distorted; a "
                "supercell is advisable (olla-dft supercell).")
        if dt_fs > 2.0:
            report.append(
                f"  WARNING: dt = {dt_fs:g} fs is large. With hydrogen or stiff "
                "bonds the total energy drifts; 0.5-1 fs is usual.")

    common = dict(
        atoms=work_atoms,
        pseudos=pseudos,
        prefix=prefix,
        pseudo_dir=pseudo_dir,
        ecutwfc=ecutwfc,
        ecutrho=ecutrho,
        insulator=opts.insulator,
        degauss=degauss,
        smearing=smearing,
        nspin=opts.nspin,
        magnetization=opts.magnetization,
        vdw=opts.vdw,
        soc=opts.soc,
        hubbard=opts.hubbard,
        hubbard_style=opts.hubbard_style,
        tot_charge=opts.tot_charge,
        dipole_correction=opts.dipole,
        hibrido=opts.hibrido,
        exx_grid=opts.exx_grid,
        exx_fraction=opts.exx_fraction,
        # La MD rompe la simetría en cuanto un átomo se mueve; si pw.x arrancó
        # con operaciones de simetría, al primer paso aborta con "some of the
        # original symmetry operations not satisfied". Por eso va forzado.
        nosym=opts.nosym or preset in ("md",),
    )

    # --- archivos por preset ---
    files = {}  # nombre -> contenido

    if preset in ("scf", "bands", "dos", "all"):
        files["scf.in"] = build_pw_input(
            calculation="scf", kcard=_kgrid_card(grid_scf), **common
        )
    if preset == "relax":
        files["relax.in"] = build_pw_input(
            calculation="relax", kcard=_kgrid_card(grid_scf), **common
        )
    if preset == "vc-relax":
        files["vc-relax.in"] = build_pw_input(
            calculation="vc-relax", kcard=_kgrid_card(grid_scf), **common
        )
    if preset == "md":
        md = dict(opts.md or {})
        files["md.in"] = build_pw_input(
            calculation="md", kcard=_kgrid_card(grid_scf), md=md, **common
        )
    if preset in ("nscf", "dos", "all"):
        files["nscf.in"] = build_pw_input(
            calculation="nscf",
            kcard=_kgrid_card(grid_nscf),
            nbnd=nbnd,
            tetrahedra=True,
            **common,
        )
        files["dos.in"] = build_dos_input(prefix)
        files["projwfc.in"] = build_projwfc_input(prefix)
        if preset == "nscf":
            files.pop("dos.in"), files.pop("projwfc.in")
    if preset in ("bands", "all"):
        kcard, labels = kpoints.kpath_card(kpath, opts.band_points)
        files["bands.in"] = build_pw_input(
            calculation="bands", kcard=kcard, nbnd=nbnd, **common
        )
        files["bands_pp.in"] = build_bandsx_input(prefix)
        label_lines = ["# index  label   kx  ky  kz  (fractional)"]
        for idx, lab, coords in labels:
            label_lines.append(
                f"{idx:4d}  {lab:8s} {coords[0]:10.6f} {coords[1]:10.6f} {coords[2]:10.6f}"
            )
        files["KPATH.txt"] = "\n".join(label_lines) + "\n"

    order = [n for n in ("scf.in", "relax.in", "vc-relax.in", "md.in", "nscf.in",
                         "dos.in", "projwfc.in", "bands.in", "bands_pp.in")
             if n in files]
    files["run.sh"] = build_run_script(order, nproc)
    files["run.py"] = build_run_python_script(order, nproc)

    # --- escritura ---
    outdir = Path(opts.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    for fname, content in files.items():
        if fname.endswith((".sh", ".py")):
            # con finales de linea POSIX y permiso de ejecucion donde exista:
            # un .sh con CRLF falla en WSL con "bad interpreter: /bin/bash^M"
            plataforma.escribir_script(outdir / fname, content)
        else:
            (outdir / fname).write_text(content, encoding="utf-8")

    report.append("")
    report.append(f"Files written to '{outdir.resolve()}':")
    report.append("  " + "  ".join(files.keys()))
    report.append("")
    report.append("Execution order (or simply ./run.sh):")
    for fname in order:
        report.append(f"  {_RUN_CMDS[fname]} -in {fname}")
    return "\n".join(report)
