# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Potenciales interatómicos aprendidos: pre-relajación y cribado.

Un MLIP (MACE, CHGNet, M3GNet...) da energías y fuerzas entre mil y diez
mil veces más barato que DFT. Aquí NO se usa para sustituir a Quantum
ESPRESSO, sino para llegar a él con el trabajo caro ya hecho:

- **pre-relajar**: entregarle a pw.x una geometría casi convergida, para
  que la relajación DFT gaste pocos pasos iónicos en vez de decenas;
- **acotar la EOS**: encontrar V0 aproximado en un segundo y colocar los
  puntos DFT centrados y en un rango estrecho, en vez de adivinarlo;
- **cribar estabilidad dinámica**: si el MLIP ya da frecuencias
  imaginarias grandes, la estructura no está relajada y la DFPT sería
  tiempo perdido — y la DFPT es lo más caro de todo Olla-DFT.

LA ADVERTENCIA QUE NO SE PUEDE OMITIR
-------------------------------------
Los modelos fundacionales están entrenados con datos PBE (o PBE+U) de
Materials Project. Si tus cálculos son LDA —o PBEsol, o con otro
funcional— el MLIP describe OTRA superficie de energía:

    Si, parámetro de red:  MACE-MP 5.464 Å | LDA 5.402 Å | exp. 5.431 Å

La diferencia no es error del modelo, es que PBE sobreestima donde LDA
subestima. Por eso:

1. una geometría relajada con MLIP es un PUNTO DE PARTIDA, nunca el
   resultado final: hay que rematar con DFT;
2. sus energías NO se pueden mezclar con las de tus cálculos DFT. Olla-DFT
   marca la procedencia y `olla-dft audit` se niega a compararlas.

DEPENDENCIA OPCIONAL
--------------------
torch y el paquete del modelo NO son requisitos de Olla-DFT. Si no están
instalados, todo lo demás funciona igual y este módulo explica qué
instalar en vez de reventar con un ImportError.
"""

import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from qekit.core.errors import ErrorDeUso

MODELOS = {
    "mace": ("mace-torch", "MACE-MP-0, trained on Materials Project (PBE)"),
    "chgnet": ("chgnet", "CHGNet, with charge states (PBE+U)"),
    "m3gnet": ("matgl", "M3GNet via matgl (PBE)"),
}


@dataclass
class MlipRun:
    modelo: str = ""
    detalle: str = ""
    funcional_entrenamiento: str = "PBE"
    atoms_inicial: object = None
    atoms_final: object = None
    pasos: int = 0
    tiempo_s: float = 0.0
    fmax_inicial: float = None
    fmax_final: float = None
    presion_inicial: float = None      # GPa
    presion_final: float = None
    desplazamiento_max: float = None      # Å
    cambio_volumen: float = None          # %
    warnings: list = field(default_factory=list)


def _falta(nombre: str, paquete: str) -> ImportError:
    return ImportError(
        f"using '{nombre}' requires installing '{paquete}', which is NOT "
        "a dependency of Olla-DFT:\n"
        f"    pip install torch {paquete}\n\n"
        "It takes a bit more than 1 GB. Everything else in Olla-DFT works without it; "
        "the MLIP only\nserves to reach the DFT calculation with the geometry "
        "almost done."
    )


def calculator(modelo: str = "mace", device: str = "cpu",
               model_size: str = "small"):
    """Devuelve un calculador de ASE para el modelo pedido."""
    modelo = (modelo or "mace").lower()
    if modelo not in MODELOS:
        raise ErrorDeUso(f"unknown model '{modelo}'. "
                         f"Options: {', '.join(MODELOS)}")
    paquete, _detalle = MODELOS[modelo]
    if modelo == "mace":
        try:
            from mace.calculators import mace_mp
        except ImportError as exc:
            raise _falta("mace", paquete) from exc
        return mace_mp(model=model_size, default_dtype="float64",
                       device=device)
    if modelo == "chgnet":
        try:
            from chgnet.model.dynamics import CHGNetCalculator
        except ImportError as exc:
            raise _falta("chgnet", paquete) from exc
        return CHGNetCalculator()
    try:
        import matgl
        from matgl.ext.ase import PESCalculator
    except ImportError as exc:
        raise _falta("m3gnet", paquete) from exc
    pot = matgl.load_model("M3GNet-MP-2021.2.8-PES")
    return PESCalculator(pot)


def relax(atoms, modelo: str = "mace", fmax: float = 0.01,
          steps: int = 300, cell: bool = True, device: str = "cpu",
          model_size: str = "small") -> MlipRun:
    """Relaja con el MLIP. `cell=True` relaja también la celda."""
    from ase.optimize import BFGS

    at = atoms.copy()
    run = MlipRun(modelo=modelo, detalle=MODELOS[modelo][1],
                  atoms_inicial=atoms.copy())
    at.calc = calculator(modelo, device=device, model_size=model_size)
    run.fmax_inicial = float(np.max(np.linalg.norm(at.get_forces(), axis=1)))
    # En un cristal simétrico las fuerzas son CERO en cualquier parámetro
    # de red: lo que mueve la celda es el esfuerzo. Reportar solo la fuerza
    # haría pensar que no pasó nada.
    try:
        run.presion_inicial = float(
            -np.trace(at.get_stress(voigt=False)) / 3.0 * 160.21766208)
    except Exception:                                  # noqa: BLE001
        pass

    objetivo = at
    if cell:
        try:
            from ase.filters import FrechetCellFilter as _Filtro
        except ImportError:                            # ASE antiguo
            from ase.constraints import UnitCellFilter as _Filtro
        objetivo = _Filtro(at)

    t0 = time.time()
    opt = BFGS(objetivo, logfile=None)
    opt.run(fmax=fmax, steps=steps)
    run.tiempo_s = time.time() - t0
    run.pasos = int(opt.get_number_of_steps())
    run.atoms_final = at
    run.fmax_final = float(np.max(np.linalg.norm(at.get_forces(), axis=1)))
    try:
        run.presion_final = float(
            -np.trace(at.get_stress(voigt=False)) / 3.0 * 160.21766208)
    except Exception:                                  # noqa: BLE001
        pass

    d = at.get_positions() - atoms.get_positions()
    run.desplazamiento_max = float(np.max(np.linalg.norm(d, axis=1)))
    v0 = atoms.get_volume()
    run.cambio_volumen = 100.0 * (at.get_volume() - v0) / v0

    if run.pasos >= steps:
        run.warnings.append(
            f"the relaxation did not reach fmax={fmax} in {steps} steps: the "
            "output geometry\nis not a minimum of the MLIP. Use it with "
            "care as a starting point.")
    if run.desplazamiento_max > 0.5:
        run.warnings.append(
            f"some atom moved {run.desplazamiento_max:.2f} Å. That is "
            "a lot: check that the\nstarting structure was what you "
            "thought, and that the MLIP describes this chemistry well.")
    return run


def volume_scan(atoms, modelo: str = "mace", span: float = 0.10,
                npoints: int = 15, device: str = "cpu",
                model_size: str = "small") -> dict:
    """Barrido E(V) con el MLIP para acotar dónde poner los puntos DFT."""
    calc = calculator(modelo, device=device, model_size=model_size)
    escalas = np.linspace(1.0 - span, 1.0 + span, npoints)
    V, E = [], []
    t0 = time.time()
    for s in escalas:
        c = atoms.copy()
        c.set_cell(atoms.cell.array * s, scale_atoms=True)
        c.calc = calc
        V.append(c.get_volume())
        E.append(c.get_potential_energy())
    V, E = np.array(V), np.array(E)

    # parábola solo para localizar el mínimo; el ajuste bueno es el de
    # Birch-Murnaghan que hace el módulo eos sobre los puntos DFT
    coef = np.polyfit(V, E, 2)
    v0 = float(-coef[1] / (2.0 * coef[0])) if coef[0] > 0 else float(V[np.argmin(E)])
    b0 = float(2.0 * coef[0] * v0 * 160.21766208) if coef[0] > 0 else None
    dentro = V.min() < v0 < V.max()
    return {
        "modelo": modelo, "V": V, "E": E, "V0": v0, "B0_aprox": b0,
        "tiempo_s": time.time() - t0,
        "escala_optima": float((v0 / atoms.get_volume()) ** (1.0 / 3.0)),
        "minimo_dentro": bool(dentro),
    }


def frequencies(atoms, modelo: str = "mace", supercell=(3, 3, 3),
                delta: float = 0.01, device: str = "cpu",
                model_size: str = "small", calc=None):
    """Frecuencias de una supercelda por diferencias finitas (cm^-1).

    Los modos de una supercelda NxNxN muestrean N^3 puntos q conmensurados
    de la zona de Brillouin: no es una dispersion, pero para la energia
    libre vibracional de la QHA es justo lo que hace falta. Devuelve
    tambien cuantas celdas primitivas contiene, para poder normalizar.
    """
    calc = calc or calculator(modelo, device=device, model_size=model_size)
    sc = atoms.repeat(tuple(int(x) for x in supercell))
    sc.calc = calc
    n = len(sc)
    masas = sc.get_masses()
    H = np.zeros((3 * n, 3 * n))
    pos0 = sc.get_positions().copy()
    for i in range(n):
        for a in range(3):
            for signo in (+1, -1):
                p = pos0.copy()
                p[i, a] += signo * delta
                sc.set_positions(p)
                H[3 * i + a] += -signo * sc.get_forces().ravel() / (2 * delta)
    sc.set_positions(pos0)
    H = 0.5 * (H + H.T)
    m = np.repeat(masas, 3)
    ev = np.linalg.eigvalsh(H / np.sqrt(np.outer(m, m)))
    freq = np.sign(ev) * np.sqrt(np.abs(ev)) * 521.4708
    celdas = int(np.prod([int(x) for x in supercell]))
    return freq, celdas


def phonon_check(atoms, modelo: str = "mace", supercell=(2, 2, 2),
                 delta: float = 0.01, device: str = "cpu",
                 model_size: str = "small") -> dict:
    """Frecuencias de TODOS los modos de la supercelda por diferencias
    finitas con el MLIP.

    Se diagonaliza la matriz dinámica completa de la supercelda (3N modos,
    con N los átomos de la supercelda), así que salen los modos de Γ de la
    celda primitiva y además los de los puntos q que la supercelda pliega
    sobre Γ. Es un CRIBADO, no un cálculo de fonones: sirve para decidir
    si vale la pena lanzar la DFPT, que puede tardar horas. Si aquí ya
    salen frecuencias imaginarias grandes, la estructura no está en un
    mínimo y la DFPT lo confirmaría después de gastar todo ese tiempo.
    """
    calc = calculator(modelo, device=device, model_size=model_size)
    sc = atoms.repeat(tuple(int(x) for x in supercell))
    sc.calc = calc
    n = len(sc)
    masas = sc.get_masses()
    t0 = time.time()

    H = np.zeros((3 * n, 3 * n))
    pos0 = sc.get_positions().copy()
    for i in range(n):
        for a in range(3):
            for signo in (+1, -1):
                p = pos0.copy()
                p[i, a] += signo * delta
                sc.set_positions(p)
                f = sc.get_forces().ravel()
                H[3 * i + a] += -signo * f / (2.0 * delta)
    sc.set_positions(pos0)
    H = 0.5 * (H + H.T)

    m = np.repeat(masas, 3)
    D = H / np.sqrt(np.outer(m, m))
    ev = np.linalg.eigvalsh(D)
    # eV/(Å²·uma) -> cm⁻¹
    CONV = 521.4708
    freq = np.sign(ev) * np.sqrt(np.abs(ev)) * CONV

    imag = freq[freq < -5.0]
    return {
        "modelo": modelo, "frecuencias": freq, "supercell": tuple(supercell),
        "n_imaginarias": int(imag.size),
        "peor_imaginaria": float(imag.min()) if imag.size else 0.0,
        "max": float(freq.max()), "tiempo_s": time.time() - t0,
        "estable": bool(imag.size == 0),
    }


MARCA = "MLIP_PROCEDENCIA.json"


def write_provenance(run: MlipRun, destino) -> str:
    """Deja junto a la estructura la marca de que la produjo un MLIP.

    Sin esto, dentro de tres meses una geometría relajada con MACE es
    indistinguible de una relajada con QE, y su energía podría acabar
    restándose contra energías DFT. `olla-dft audit` busca este archivo.
    """
    import json

    d = Path(destino)
    carpeta = d.parent if d.suffix else d
    carpeta.mkdir(parents=True, exist_ok=True)
    doc = {
        "origen": f"mlip:{run.modelo}",
        "modelo": run.detalle,
        "funcional_entrenamiento": run.funcional_entrenamiento,
        "archivo": d.name if d.suffix else None,
        "pasos": run.pasos,
        "fmax_final_eV_A": run.fmax_final,
        "presion_final_GPa": run.presion_final,
        "cambio_volumen_pct": run.cambio_volumen,
        "aviso": ("Geometry produced by a machine-learned potential, not by "
                  "DFT. Its energy is NOT comparable with Quantum "
                  "ESPRESSO energies: different energy surface."),
    }
    f = carpeta / MARCA
    f.write_text(json.dumps(doc, ensure_ascii=False, indent=2) + "\n")
    return str(f)


def read_provenance(carpeta):
    """Lee la marca de procedencia MLIP de una carpeta, si existe."""
    import json

    f = Path(carpeta)
    f = (f if f.is_dir() else f.parent) / MARCA
    if not f.exists():
        return None
    try:
        return json.loads(f.read_text())
    except Exception:                                  # noqa: BLE001
        return None


# ----------------------------------------------------------------------
def report_relax(run: MlipRun, funcional_dft: str = None) -> str:
    lines = ["--- Pre-relaxation with machine-learned potential ---",
             f"Model: {run.detalle}",
             f"{run.pasos} steps in {run.tiempo_s:.2f} s",
             f"Maximum force: {run.fmax_inicial:.4f} -> "
             f"{run.fmax_final:.4f} eV/Å",]
    if run.presion_inicial is not None:
        lines.append(
            f"Pressure: {run.presion_inicial:+.2f} -> "
            f"{run.presion_final:+.2f} GPa   (in a symmetric crystal the "
            "forces are\n  zero by construction; what relaxes the cell "
            "is the stress)")
    lines += [
             f"Maximum displacement: {run.desplazamiento_max:.4f} Å  |  "
             f"volume change: {run.cambio_volumen:+.2f} %"]
    for w in run.warnings:
        lines.append(f"\nWARNING: {w}")
    lines += ["",
              "THIS IS NOT THE FINAL RESULT. The model is trained on "
              "PBE data from\nMaterials Project; if your calculation uses another "
              "functional, it describes a different\nenergy surface. "
              "Use this geometry as the STARTING POINT of a DFT relax,\n"
              "and do not mix its energies with those of QE."]
    if funcional_dft and funcional_dft.upper() not in ("PBE", "PBESOL"):
        lines.append(
            f"\nYour calculation uses {funcional_dft}: the difference from PBE will be "
            "systematic.\nIn silicon, for example, MACE gives 5.464 Å and "
            "LDA 5.402 Å.")
    return "\n".join(lines)


def report_scan(d: dict, atoms) -> str:
    lines = ["--- Volume scan with machine-learned potential ---",
             f"Model: {MODELOS[d['modelo']][1]}",
             f"{len(d['V'])} points in {d['tiempo_s']:.2f} s",
             f"Approximate V0: {d['V0']:.3f} Å³  "
             f"(scale {d['escala_optima']:.4f} on the given cell)"]
    if d["B0_aprox"]:
        lines.append(f"Approximate B0: {d['B0_aprox']:.1f} GPa")
    if not d["minimo_dentro"]:
        lines.append(
            "\nWARNING: the minimum fell OUTSIDE the scanned range. Widen --span; "
            "the V0 above\nis an extrapolation and is not suitable to centre "
            "the DFT points.")
    else:
        lines.append(
            "\nUse the optimal scale to centre the DFT EOS:\n"
            f"    olla-dft eos structure.cif --scale {d['escala_optima']:.4f}"
            " --span 0.04\n"
            "With the minimum already located, a narrow range gives a "
            "better fit with fewer\npoints — which is where the real saving is.")
    lines.append(
        "\nThe B0 here comes from a parabola and from the model's functional: "
        "it serves to know\nthe order of magnitude, not to report it.")
    return "\n".join(lines)


def report_phonon(d: dict) -> str:
    lines = ["--- Dynamical stability screening (MLIP) ---",
             f"Model: {MODELOS[d['modelo']][1]}",
             f"Supercell {d['supercell'][0]}x{d['supercell'][1]}x"
             f"{d['supercell'][2]}  |  {d['tiempo_s']:.1f} s",
             f"Maximum frequency: {d['max']:.1f} cm⁻¹"]
    if d["estable"]:
        lines += ["",
                  "No imaginary frequencies: the structure appears to be "
                  "at a minimum.\nIt is worth launching the DFPT."]
    else:
        lines += ["",
                  f"{d['n_imaginarias']} imaginary frequencies, the worst "
                  f"at {d['peor_imaginaria']:.1f} cm⁻¹.",
                  "The structure is NOT at a minimum of the model. Before "
                  "spending hours on DFPT,\nrelax it better (with --pre-ml or with "
                  "a tighter DFT vc-relax)."]
    lines.append(
        "\nThis is a SCREENING, not a phonon calculation: the model is not your "
        "functional and\nfinite differences on a small supercell "
        "do not reproduce the\ndispersion. It serves to decide whether to launch the "
        "DFPT, not to replace it.")
    return "\n".join(lines)
