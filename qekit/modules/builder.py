# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Constructores de estructuras: superficies y defectos puntuales.

Dos cosas que hacen falta constantemente y que nadie quiere hacer a mano:

**Superficies.** Cortar un cristal por unos índices de Miller (hkl), apilar
N capas y añadir vacío. Olla-DFT usa el constructor de ASE y encima le pone lo
que un cálculo necesita de verdad: congelar las capas del fondo para que
simulen el volumen, avisar si el vacío es insuficiente, y contar si la losa
es simétrica (las dos caras iguales) o polar (caras distintas, con dipolo).

**Defectos puntuales.** Vacancia, sustitución e intersticial dentro de una
supercelda, con la fórmula de energía de formación lista para llenar.

Detalles que importan:
- una losa POLAR genera un dipolo artificial a través del vacío por las
  condiciones periódicas; el reporte lo detecta y recuerda `dipfield`;
- el vacío se mide entre superficies ATÓMICAS, no entre bordes de celda,
  que es donde se cuela el error de "puse 15 Å" cuando en realidad hay 9;
- para un defecto, la energía de formación depende del potencial químico
  del reservorio, y no hay un valor "correcto" universal: Olla-DFT escribe la
  fórmula con sus términos y deja explícito lo que hay que elegir.
"""

from dataclasses import dataclass, field

import numpy as np
from ase import Atoms

from qekit.core import structure as struct_mod
from qekit.core.errors import ErrorDeUso


@dataclass
class SlabInfo:
    atoms: Atoms = None
    miller: tuple = None
    layers: int = 0
    vacuum: float = 0.0
    vacuum_real: float = 0.0        # hueco entre superficies atómicas
    thickness: float = 0.0
    simetrica: bool = False
    polar: bool = False
    fijados: int = 0
    planos: int = 0                 # planos atómicos distintos en z
    planos_fijos: int = 0
    warnings: list = field(default_factory=list)


def surface(atoms: Atoms, miller=(1, 0, 0), layers: int = 4,
            vacuum: float = 15.0, fix_layers: int = 0,
            tol: float = 0.3) -> SlabInfo:
    """Corta una superficie (hkl) con vacío, sobre la celda convencional.

    `fix_layers` congela ese número de capas atómicas del fondo (se marcan
    en el input de pw.x para que simulen el volumen y no se relajen).

    Hubo aquí un parámetro `orthogonal` que no hacía NADA: aceptaba el
    argumento, entraba en un bloque vacío y devolvía la misma losa. Nadie lo
    usaba, así que se quitó en vez de implementarlo. Una bandera que se
    acepta y se ignora es peor que no tenerla: el usuario cree que pidió
    algo.
    """
    from ase.build import surface as ase_surface

    base = struct_mod.conventional(atoms)
    miller = tuple(int(m) for m in miller)
    if all(m == 0 for m in miller):
        raise ErrorDeUso("the Miller indices cannot be (0,0,0)")

    slab = ase_surface(base, miller, layers, vacuum=vacuum / 2.0,
                       periodic=True)
    slab.center(vacuum=vacuum / 2.0, axis=2)

    z = slab.positions[:, 2]
    grosor = float(z.max() - z.min())
    c = float(slab.cell.array[2, 2])
    vac_real = c - grosor

    info = SlabInfo(atoms=slab, miller=miller, layers=layers, vacuum=vacuum,
                    vacuum_real=vac_real, thickness=grosor)
    # el corte se hace sobre la celda CONVENCIONAL, que es la referencia de
    # los índices de Miller. El precio es que la celda superficial puede
    # salir como un múltiplo de la mínima: conviene decirlo, porque el costo
    # del cálculo va con el número de átomos por capa.
    por_plano = len(slab) / max(len(_planos_z(slab, tol)), 1)
    if por_plano > 1.5:
        info.warnings.append(
            f"the surface cell has {por_plano:.0f} atoms per plane. "
            "The cut is made on the conventional cell (which is the "
            "reference for the hkl indices) and that can give a multiple of "
            "the minimal surface cell. If you only need the clean surface, "
            "a smaller cell is cheaper.")

    # ¿es simétrica? se compara el perfil de z con su reflejo
    zc = np.sort(z - z.mean())
    info.simetrica = bool(np.allclose(zc, -zc[::-1], atol=tol))
    # polar: composición de la capa superior distinta de la inferior
    sup = _capa(slab, z.max(), tol)
    inf = _capa(slab, z.min(), tol)
    info.polar = sorted(sup) != sorted(inf)

    info.planos = len(_planos_z(slab, tol))
    if fix_layers > 0:
        info.fijados = _fijar_capas(slab, fix_layers, tol)
        info.planos_fijos = min(fix_layers, info.planos)
        if fix_layers >= info.planos:
            info.warnings.append(
                f"you asked to freeze {fix_layers} planes and the slab only has "
                f"{info.planos}: it would be entirely fixed and there would be "
                "no surface relaxation, which is exactly what one wants to compute.")

    if vac_real < 10.0:
        info.warnings.append(
            f"the REAL vacuum between atomic surfaces is {vac_real:.1f} Å "
            f"(you asked for {vacuum:.1f} of cell). Below ~10 Å the two "
            "faces see each other and the surface energy and work function "
            "come out wrong.")
    if info.polar:
        info.warnings.append(
            "the slab is POLAR (the two faces are not equivalent): the "
            "periodic boundary conditions create an artificial dipole across "
            "the vacuum. Add 'dipfield = .true.' and 'edir = 3' to the input, "
            "or cut a symmetric slab.")
    if layers < 4:
        info.warnings.append(
            f"{layers} layers is too few: the centre of the slab should "
            "resemble the bulk, and with so few layers it does not.")
    return info


def _capa(slab: Atoms, z0: float, tol: float) -> list:
    sel = np.abs(slab.positions[:, 2] - z0) < tol
    return [s for s, m in zip(slab.get_chemical_symbols(), sel) if m]


def _planos_z(slab: Atoms, tol: float) -> list:
    """Alturas z distintas ocupadas por átomos (los planos atómicos).

    OJO: un "plano atómico" no es lo mismo que una "capa" del constructor
    de ASE. En el diamante (100), por ejemplo, cada capa de ASE contiene
    varios planos de z distintos. Olla-DFT cuenta planos, que es lo que
    importa para decidir qué congelar, y reporta cuántos hay.
    """
    z = slab.positions[:, 2]
    niveles = []
    for zi in np.sort(z):
        if not niveles or abs(zi - niveles[-1]) > tol:
            niveles.append(float(zi))
    return niveles


#: Formato de estructura que conserva los átomos fijos al exportar. El CIF
#: no tiene dónde guardarlos; POSCAR sí ("Selective dynamics") y ASE los
#: lee de vuelta como restricción FixAtoms.
FORMATO_CON_FIJOS = "POSCAR (or .vasp)"


def _fijar_capas(slab: Atoms, n: int, tol: float) -> int:
    """Marca como fijos los átomos de los n planos atómicos inferiores.

    Se guarda de dos formas: en slab.arrays['qekit_fijo'] (1 = fijo) y como
    restricción FixAtoms de ASE. `inputgen.build_pw_input` lee cualquiera
    de las dos y pone '0 0 0' en la tercera columna de ATOMIC_POSITIONS.
    La restricción FixAtoms es la que sobrevive a un POSCAR; el CIF pierde
    las dos.
    """
    from ase.constraints import FixAtoms

    z = slab.positions[:, 2]
    niveles = _planos_z(slab, tol)
    corte = niveles[min(n, len(niveles)) - 1] + tol
    fijo = (z <= corte).astype(int)
    slab.set_array("qekit_fijo", fijo)
    slab.set_constraint(FixAtoms(indices=np.flatnonzero(fijo)))
    return int(fijo.sum())


# ----------------------------------------------------------------------
@dataclass
class DefectInfo:
    atoms: Atoms = None
    kind: str = ""
    site: int = None
    especie_ida: str = ""
    especie_nueva: str = ""
    supercell: tuple = None
    n_perfecto: int = 0
    warnings: list = field(default_factory=list)


def defect(atoms: Atoms, kind: str = "vacancy", site: int = 0,
           new_element: str = None, supercell=(2, 2, 2),
           position=None) -> tuple:
    """Crea un defecto puntual en una supercelda.

    kind: 'vacancy' (quita el átomo `site`), 'substitution' (lo cambia por
    `new_element`) o 'interstitial' (mete `new_element` en `position`, en
    coordenadas fraccionarias de la supercelda).

    Devuelve (perfecto, con_defecto) para que la energía de formación se
    calcule con los dos con EXACTAMENTE la misma celda y malla.
    """
    base = struct_mod.primitive(atoms)
    n1, n2, n3 = (int(x) for x in supercell)
    perfecto = base.repeat((n1, n2, n3))
    d = perfecto.copy()
    info = DefectInfo(kind=kind, supercell=(n1, n2, n3),
                      n_perfecto=len(perfecto))

    if kind == "vacancy":
        if not 0 <= site < len(d):
            raise ErrorDeUso(f"site index {site} out of range "
                             f"(0..{len(d)-1})")
        info.especie_ida = d.get_chemical_symbols()[site]
        info.site = site
        del d[site]
    elif kind == "substitution":
        if not new_element:
            raise ErrorDeUso("the substitution needs --new-element: specify "
                             "which species goes in, e.g. --new-element P")
        if not 0 <= site < len(d):
            raise ErrorDeUso(f"site index {site} out of range "
                             f"(0..{len(d)-1})")
        info.especie_ida = d.get_chemical_symbols()[site]
        info.especie_nueva = new_element
        info.site = site
        d[site].symbol = new_element
    elif kind == "interstitial":
        if not new_element:
            raise ErrorDeUso("the interstitial needs --new-element: specify "
                             "which species is inserted, e.g. --new-element H")
        if position is None:
            raise ErrorDeUso("the interstitial needs --position x,y,z "
                             "(fractional coordinates of the supercell)")
        pos_cart = np.asarray(position, dtype=float) @ d.cell.array
        d.append(Atoms(new_element, positions=[pos_cart])[0])
        info.especie_nueva = new_element
        # ¿quedó demasiado cerca de un vecino?
        dist = d.get_distances(len(d) - 1, range(len(d) - 1), mic=True)
        if dist.min() < 1.0:
            info.warnings.append(
                f"the interstitial ended up {dist.min():.2f} Å from its nearest "
                "neighbour: check the position, that will not converge.")
    else:
        raise ErrorDeUso("kind must be vacancy, substitution or interstitial")

    info.atoms = d
    lado = min(np.linalg.norm(perfecto.cell.array, axis=1))
    if lado < 10.0:
        info.warnings.append(
            f"the supercell measures {lado:.1f} Å along its shortest side: the "
            "defect sees its periodic images. For formation energies "
            "≥ 10-12 Å is advisable.")
    return perfecto, info


def formation_energy_text(info: DefectInfo) -> str:
    """La fórmula con sus términos, para no aplicarla a ciegas."""
    if info.kind == "vacancy":
        term = f"+ mu({info.especie_ida})"
        quita = f"one {info.especie_ida} was removed"
    elif info.kind == "substitution":
        term = f"+ mu({info.especie_ida}) - mu({info.especie_nueva})"
        quita = f"{info.especie_ida} -> {info.especie_nueva}"
    else:
        term = f"- mu({info.especie_nueva})"
        quita = f"one {info.especie_nueva} was added"
    return "\n".join([
        "Formation energy:",
        f"  E_f = E(defect) - E(perfect) {term}  [+ q(E_F + E_v) + E_corr]",
        f"  ({quita})",
        "",
        "  mu = chemical potential of the reservoir. There is NO universal value:",
        "  it depends on the synthesis conditions (rich or poor in each",
        "  species) and bounds E_f between two limits instead of fixing a number.",
        "  The terms in brackets only apply to CHARGED defects.",
    ])


def report_slab(info: SlabInfo) -> str:
    a = info.atoms
    lines = ["--- Surface ---",
             f"Miller indices: ({info.miller[0]}{info.miller[1]}"
             f"{info.miller[2]})  |  {info.layers} layers  |  "
             f"{len(a)} atoms",
             f"Formula: {a.get_chemical_formula()}",
             f"Slab thickness: {info.thickness:.2f} Å",
             f"Real vacuum between surfaces: {info.vacuum_real:.2f} Å",
             f"Atomic planes along z: {info.planos}",
             f"Symmetric slab: {'yes' if info.simetrica else 'no'}  |  "
             f"polar: {'yes' if info.polar else 'no'}"]
    if info.fijados:
        lines.append(f"Frozen: {info.planos_fijos} bottom planes "
                     f"({info.fijados} atoms of {len(a)})")
    for w in info.warnings:
        lines.append(f"\nWARNING: {w}")
    return "\n".join(lines)


def report_defect(info: DefectInfo) -> str:
    a = info.atoms
    lines = ["--- Point defect ---",
             f"Type: {info.kind}  |  supercell "
             f"{info.supercell[0]}x{info.supercell[1]}x{info.supercell[2]}",
             f"Perfect: {info.n_perfecto} atoms  ->  with defect: "
             f"{len(a)} atoms ({a.get_chemical_formula()})"]
    if info.site is not None:
        lines.append(f"Affected site: index {info.site} "
                     f"({info.especie_ida})")
    if info.especie_nueva:
        lines.append(f"Introduced species: {info.especie_nueva}")
    lines += ["", formation_energy_text(info)]
    for w in info.warnings:
        lines.append(f"\nWARNING: {w}")
    return "\n".join(lines)
