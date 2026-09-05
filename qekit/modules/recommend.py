# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Recomendaciones a partir de TU propio historial de cálculos.

La base de `olla-dft db` acumula, sin que cueste nada extra, lo que hace falta
para no repetir el mismo tanteo: qué cutoffs convergieron, qué mezcla
funcionó, cuánto tardó cada cosa y en qué sistemas.

POR QUÉ ESTO NO ES UNA RED NEURONAL
-----------------------------------
Con unas decenas de cálculos —que es lo que va a haber— un modelo aprendido
sobreajusta y no se puede auditar. Aquí se usa lo que sí funciona con pocos
datos: buscar los cálculos PARECIDOS al que quieres hacer (mismos
elementos, sistema parecido) y mirar qué les funcionó, diciendo siempre
cuántos casos respaldan cada número.

Una recomendación con un solo caso detrás se marca como tal. Es la
diferencia entre "esto suele funcionar" y "esto funcionó una vez".

Y una cosa que NO hace: inventar cutoffs. Si no hay historial del elemento,
lo dice y remite a los cutoffs que declara el propio pseudopotencial, que
es un dato y no una predicción.
"""

from dataclasses import dataclass

import numpy as np


@dataclass
class Sugerencia:
    campo: str = ""
    valor: object = None
    n_casos: int = 0
    rango: tuple = None
    razon: str = ""
    confianza: str = "baja"      # "baja" | "media" | "alta"


def _confianza(n: int) -> str:
    if n >= 8:
        return "alta"
    if n >= 3:
        return "media"
    return "baja"


def similares(filas: list, elementos, natoms: int = None,
              tol_natoms: float = 2.0) -> list:
    """Cálculos del historial parecidos al que se quiere hacer.

    "Parecido" = comparte al menos un elemento; se prioriza compartir
    todos. No se usa una métrica sofisticada a propósito: con pocos datos,
    cualquier cosa más elaborada da una falsa sensación de precisión.
    """
    els = set(elementos)
    out = []
    for f in filas:
        if not f.get("convergido"):
            continue
        formula = f.get("formula") or ""
        f_els = set(_elementos_de(formula))
        if not (f_els & els):
            continue
        puntaje = len(f_els & els) / max(len(f_els | els), 1)
        if natoms and f.get("natoms"):
            razon = f["natoms"] / natoms
            if razon > tol_natoms or razon < 1.0 / tol_natoms:
                puntaje *= 0.5
        out.append((puntaje, f))
    return [f for _p, f in sorted(out, key=lambda t: -t[0])]


def _elementos_de(formula: str) -> list:
    import re
    return re.findall(r"[A-Z][a-z]?", formula or "")


def sugerir(filas: list, elementos, natoms: int = None,
            es_losa: bool = False) -> list:
    """Sugerencias de parámetros para un cálculo nuevo."""
    sug = []
    vecinos = similares(filas, elementos, natoms)
    if not vecinos:
        sug.append(Sugerencia(
            campo="(sin historial)", n_casos=0, confianza="baja",
            razon="There are no previous calculations with these elements. Use the "
                  "cutoffs declared\nby the pseudopotential itself (Olla-DFT "
                  "reads them from the UPF) or the SSSP table: that is a\nmeasured "
                  "datum, not a prediction, and it always beats an "
                  "extrapolation."))
        return sug

    ecuts = [f["ecutwfc"] for f in vecinos if f.get("ecutwfc")]
    if ecuts:
        v = float(np.max(ecuts))
        sug.append(Sugerencia(
            campo="ecutwfc", valor=v, n_casos=len(ecuts),
            rango=(float(np.min(ecuts)), float(np.max(ecuts))),
            confianza=_confianza(len(ecuts)),
            razon=f"the MAXIMUM of {len(ecuts)} converged calculations with "
                  f"these elements (range {min(ecuts):.0f}–{max(ecuts):.0f} "
                  "Ry). The maximum is taken, not the mean: a low cutoff that "
                  "worked in one\nsystem guarantees nothing in another."))

    duales = [(f["ecutrho"] / f["ecutwfc"]) for f in vecinos
              if f.get("ecutrho") and f.get("ecutwfc")]
    if duales:
        sug.append(Sugerencia(
            campo="dual (ecutrho/ecutwfc)", valor=float(np.max(duales)),
            n_casos=len(duales), confianza=_confianza(len(duales)),
            razon="the dual used by the previous calculations; it depends on the "
                  "pseudopotential type\n(4 for norm-conserving, 8-12 for "
                  "ultrasoft and PAW)."))

    dens = [f["kdensity"] for f in vecinos if f.get("kdensity")]
    if dens:
        sug.append(Sugerencia(
            campo="k-point density (points/Å⁻³)", valor=float(np.median(dens)),
            n_casos=len(dens), rango=(float(np.min(dens)),
                                      float(np.max(dens))),
            confianza=_confianza(len(dens)),
            razon="median of the converged calculations. The density, not the "
                  "number of points,\nis what is comparable between cells of "
                  "different size."))

    pasos = [f["n_scf"] for f in vecinos if f.get("n_scf")]
    if pasos and np.median(pasos) > 40:
        sug.append(Sugerencia(
            campo="electron_maxstep", valor=300, n_casos=len(pasos),
            confianza=_confianza(len(pasos)),
            razon=f"your calculations with these elements needed a "
                  f"median of {np.median(pasos):.0f} SCF steps:\nthe "
                  "default maximum falls short."))

    if es_losa:
        sug.append(Sugerencia(
            campo="mixing_beta", valor=0.3, n_casos=0, confianza="baja",
            razon="it is a slab with vacuum: these give the most charge "
                  "sloshing. Starting\nwith a low mixing_beta and "
                  "mixing_mode='local-TF' saves retries. This does not come "
                  "from\nyour history, it is a general rule."))
    return sug


def report(sug: list, elementos, n_historial: int = 0) -> str:
    lines = ["--- Suggestions from your history ---",
             f"Elements: {', '.join(elementos)}  |  "
             f"calculations in the database: {n_historial}", ""]
    if not sug or sug[0].campo == "(sin historial)":
        lines.append(sug[0].razon if sug else "No data.")
        return "\n".join(lines)

    for s in sug:
        val = s.valor
        if isinstance(val, float):
            val = f"{val:.4g}"
        # la confianza 'baja' cubre 1 y 2 casos: no se dice "un solo caso"
        # cuando hay dos
        marca = {"alta": "", "media": "  (few cases)",
                 "baja": ("  (A SINGLE CASE: take it as a hint)"
                          if s.n_casos == 1 else
                          f"  (ONLY {s.n_casos} CASES: take it as a hint)")
                 }[s.confianza]
        if s.n_casos == 0:
            marca = "  (general rule, not from your history)"
        lines.append(f"  {s.campo}: {val}"
                     f"   [{s.n_casos} case{'s' if s.n_casos != 1 else ''}]"
                     f"{marca}")
        for l in s.razon.splitlines():
            lines.append(f"      {l}")
        lines.append("")
    lines.append("These suggestions come from what ALREADY worked for you, not from a "
                 "trained model.\nThey do not replace a convergence "
                 "test: 'olla-dft converge' is still\nthe way to "
                 "know for sure for a new system.")
    return "\n".join(lines)
