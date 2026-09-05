# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Cuánto va a tardar esto, antes de lanzarlo.

La pregunta que se hace todo el mundo antes de darle a --run es si el barrido
cabe en la tarde. Hoy se responde a ojo, y a ojo se falla por un factor diez.

La idea aquí no es predecir el tiempo de pw.x en abstracto —eso depende de la
máquina, de la compilación, de si hay otra cosa corriendo— sino calibrar un
modelo físico sencillo con los cálculos que TÚ ya has hecho, que están en la
base de `olla-dft db`. El modelo pone la forma y el historial pone la escala.

La forma sale de cómo escala un cálculo de ondas planas:

    t ≈ C · n_k · n_espín · N_PW · n_bandas · n_scf

  N_PW    número de ondas planas, = V·ecut^(3/2)/(6π²) en unidades atómicas.
          Se calcula de verdad, no se aproxima con una potencia del cutoff:
          así el mismo modelo vale para una celda pequeña con cutoff alto y
          para una supercelda con cutoff bajo.
  n_k     puntos k IRREDUCIBLES, contados con spglib. Usar la malla entera
          erraría por un factor de hasta 48 en un cristal cúbico.
  n_scf   iteraciones, que no se saben de antemano: se toma la mediana de tu
          propio historial y se dice cuál se usó.

Lo que se devuelve siempre lleva su dispersión. Un número solo, sin decir que
tus cálculos anteriores se reparten en un factor dos alrededor del modelo,
sería más cómodo y menos cierto.
"""

import math
import re
import sqlite3
import statistics
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np


BOHR = 0.529177210903
# 1 Å³ en bohr³
A3_BOHR3 = 1.0 / BOHR ** 3

# iteraciones scf típicas si no hay historial del que sacarlas
N_SCF_DEFECTO = 14
# pasos iónicos típicos de una relajación, si no hay historial
PASOS_IONICOS = {"relax": 8.0, "vc-relax": 12.0, "md": 1.0, "scf": 1.0,
                 "nscf": 1.0, "bands": 1.0}


def n_ondas_planas(volumen_A3: float, ecutwfc_Ry: float) -> float:
    """N_PW = V·k_max³/(6π²) con k_max = √ecut en unidades atómicas de Ry."""
    if not volumen_A3 or not ecutwfc_Ry or ecutwfc_Ry <= 0:
        return 0.0
    v_bohr = float(volumen_A3) * A3_BOHR3
    return v_bohr * float(ecutwfc_Ry) ** 1.5 / (6.0 * math.pi ** 2)


def k_irreducibles(atoms, malla, shift=(0, 0, 0)) -> int:
    """Puntos k irreducibles de la malla, con la simetría de la estructura.

    Sin esto el modelo erraría por el factor de reducción, que va de 1 en una
    celda sin simetría a 48 en una cúbica: exactamente el rango en el que la
    estimación deja de servir para nada.
    """
    try:
        import spglib
        celda = (np.asarray(atoms.cell.array), atoms.get_scaled_positions(),
                 atoms.get_atomic_numbers())
        mapa, _ = spglib.get_ir_reciprocal_mesh(
            list(int(x) for x in malla), celda, is_shift=list(shift))
        return int(len(np.unique(mapa)))
    except Exception:                                       # noqa: BLE001
        return int(np.prod([int(x) for x in malla]))


# ----------------------------------------------------------------------
# Descriptores de un cálculo, leídos de su propio input
# ----------------------------------------------------------------------
def _valor(texto: str, clave: str):
    m = re.search(rf"^\s*{clave}\s*=\s*([^\s,!]+)", texto,
                  re.IGNORECASE | re.MULTILINE)
    if not m:
        return None
    v = m.group(1).strip().strip("'\"")
    if v.lower() in (".true.", ".false."):
        return v.lower() == ".true."
    try:
        return float(v)
    except ValueError:
        return v


_RE_NK = re.compile(r"number of k points\s*=\s*(\d+)")


def nk_de_salida(carpeta) -> int:
    """Puntos k que pw.x dijo usar, si esa carpeta ya corrió alguna vez.

    Vale más que cualquier estimación: pw.x y spglib no siempre encuentran
    las mismas operaciones de simetría. En silicio en celda primitiva,
    spglib ve 48 operaciones y una malla 9x9x9 le da 35 puntos
    irreducibles, mientras que pw.x encuentra 12 operaciones y usa 85. Un
    factor 2.4 en el coste previsto sale de ahí, no de la física.
    """
    for nombre in ("pw.out", "scf.out", "espresso.out"):
        f = Path(carpeta) / nombre
        if not f.exists():
            continue
        try:
            m = _RE_NK.search(f.read_text(errors="ignore"))
        except OSError:
            continue
        if m:
            return int(m.group(1))
    return 0


def descriptores_de_input(ruta) -> dict:
    """Lee de un pw.in lo que hace falta para estimar su coste."""
    from qekit.core import structure

    ruta = Path(ruta)
    texto = ruta.read_text(errors="ignore")
    d = {"ruta": str(ruta)}
    d["calculation"] = _valor(texto, "calculation") or "scf"
    d["ecutwfc"] = _valor(texto, "ecutwfc") or 0.0
    d["nat"] = int(_valor(texto, "nat") or 0)
    d["nspin"] = int(_valor(texto, "nspin") or 1)
    if _valor(texto, "noncolin") is True:
        d["nspin"] = 4          # no colineal: espinores, cuatro componentes
    nbnd = _valor(texto, "nbnd")

    atoms = None
    try:
        atoms = structure.load(str(ruta))
    except Exception:                                       # noqa: BLE001
        atoms = None
    d["volumen"] = (float(abs(np.linalg.det(atoms.cell.array)))
                    if atoms is not None else 0.0)
    d["formula"] = (atoms.get_chemical_formula() if atoms is not None else "")
    if not d["nat"] and atoms is not None:
        d["nat"] = len(atoms)

    # bandas: las del input si están, y si no una estimación por electrones
    if nbnd:
        d["nbnd"] = float(nbnd)
    else:
        d["nbnd"] = max(4.0, d["nat"] * 2.0)

    # puntos k
    m = re.search(r"^\s*K_POINTS\s*(\{?\s*\w*\s*\}?)\s*$(.*)", texto,
                  re.IGNORECASE | re.MULTILINE | re.DOTALL)
    d["nk"] = 1
    if m:
        modo = (m.group(1) or "").strip("{} \t").lower()
        resto = m.group(2).strip().splitlines()
        if modo.startswith("automatic") and resto:
            nums = [int(float(x)) for x in resto[0].split()[:6]]
            malla = nums[:3] or [1, 1, 1]
            shift = nums[3:6] if len(nums) >= 6 else [0, 0, 0]
            d["malla"] = tuple(malla)
            d["nk"] = (k_irreducibles(atoms, malla, shift)
                       if atoms is not None else int(np.prod(malla)))
        elif modo.startswith("gamma"):
            d["nk"] = 1
        elif resto:
            try:
                d["nk"] = max(1, int(resto[0].split()[0]))
            except (ValueError, IndexError):
                d["nk"] = 1
    d["nk_fuente"] = "symmetry estimate (spglib)"
    real = nk_de_salida(ruta.parent)
    if real:
        d["nk"] = real
        d["nk_fuente"] = "the one pw.x used"
    d["npw"] = n_ondas_planas(d["volumen"], d["ecutwfc"])
    return d


def trabajo(d: dict) -> tuple:
    """Los dos términos del coste de una iteración scf.

    Un paso de pw.x tiene dos partes que escalan distinto:

      w1 = n_k · espín · N_PW · n_bandas          las FFT, una por banda
      w2 = n_k · espín · N_PW · n_bandas²         ortogonalizar y rotar el
                                                  subespacio, que va con el
                                                  cuadrado de las bandas

    En una celda pequeña manda w1 y en una supercelda manda w2. Con un solo
    término el modelo se queda corto justo donde importa: predijo 65 s para
    un cálculo de 53 átomos que tardó 265. Se devuelven los dos y el ajuste
    decide el peso de cada uno con el historial.
    """
    espin = {1: 1.0, 2: 2.0, 4: 4.0}.get(int(d.get("nspin", 1) or 1), 1.0)
    base = max(1e-12, d.get("nk", 1) * espin * d.get("npw", 0.0))
    nb = max(1.0, float(d.get("nbnd", 1.0)))
    return base * nb, base * nb * nb


def iteraciones(d: dict, modelo=None) -> float:
    """Iteraciones scf totales previstas, contando los pasos iónicos."""
    calc = d.get("calculation", "scf")
    if modelo is not None and modelo.calibrado:
        n_scf = modelo.n_scf_mediana
        iones = modelo.pasos_ionicos.get(calc, PASOS_IONICOS.get(calc, 1.0))
    else:
        n_scf = N_SCF_DEFECTO
        iones = PASOS_IONICOS.get(calc, 1.0)
    return max(1.0, n_scf) * max(1.0, iones)


# ----------------------------------------------------------------------
# Calibración con el historial
# ----------------------------------------------------------------------
@dataclass
class Modelo:
    C1: float = None                # segundos por FFT y banda
    C2: float = None                # segundos por ortogonalización
    t0: float = 0.0                 # segundos fijos de arranque por cálculo
    n: int = 0                      # cálculos usados
    dispersion: float = None        # factor multiplicativo (1.6 = ±60 %)
    n_scf_mediana: float = N_SCF_DEFECTO
    pasos_ionicos: dict = field(default_factory=dict)
    rango: float = 1.0              # cuánto varía el trabajo en el historial
    sesgo: float = None             # sesgo medido fuera de muestra
    n_sistemas: int = 0             # combinaciones distintas en el historial
    validado: bool = False          # ¿la dispersión sale de validación cruzada?
    fuente: str = ""

    @property
    def calibrado(self) -> bool:
        return self.C1 is not None and self.n > 0

    @property
    def extrapola_bien(self) -> bool:
        """¿Da el historial para predecir un sistema distinto?

        Si todos los cálculos indexados son casi iguales, la constante sale
        bien ajustada y aun así la predicción de un sistema tres veces más
        grande puede fallar por un factor dos: no hay nada en los datos que
        fije la pendiente. Se midió: con diecisiete cálculos de silicio
        idénticos, predecir el mismo silicio con otro cutoff se fue un 44 %.
        """
        return self.n >= 8 and self.rango >= 5.0


def historial(db_path="olla-dft.db") -> list:
    """Cálculos con tiempo de pared y descriptores suficientes."""
    p = Path(db_path)
    if not p.exists():
        return []
    con = sqlite3.connect(str(p))
    try:
        cur = con.execute(
            "SELECT natoms, ecutwfc, kgrid, nspin, volumen_A3, n_scf, "
            "nk, nbnd, n_bfgs, wall_s, calculation FROM calculos "
            "WHERE wall_s IS NOT NULL AND wall_s > 0 AND ecutwfc > 0 "
            "AND volumen_A3 > 0")
        filas = [dict(zip([c[0] for c in cur.description], f))
                 for f in cur.fetchall()]
    except sqlite3.Error:
        return []
    finally:
        con.close()
    return filas


def _fila_a_descriptor(f: dict) -> dict:
    nat = int(f.get("natoms") or 0)
    nk = f.get("nk")
    if not nk:
        # base antigua sin la columna: se usa la malla entera, sabiendo que
        # sobreestima. Es la aproximación menos mala con lo que hay guardado.
        malla = (1, 1, 1)
        if f.get("kgrid"):
            nums = re.findall(r"\d+", str(f["kgrid"]))
            if len(nums) >= 3:
                malla = tuple(int(x) for x in nums[:3])
        nk = int(np.prod(malla))
    return {"nat": nat,
            "nbnd": float(f.get("nbnd") or max(4.0, nat * 2.0)),
            "nspin": int(f.get("nspin") or 1),
            "npw": n_ondas_planas(f.get("volumen_A3"), f.get("ecutwfc")),
            "nk": max(1, int(nk)),
            "calculation": f.get("calculation") or "scf",
            "n_bfgs": int(f.get("n_bfgs") or 1)}


def nk_de_historial(db_path, formula: str, malla) -> int:
    """nk real de un cálculo anterior con la misma fórmula y la misma malla.

    Es mejor fallback que spglib porque es literalmente lo que pw.x hizo con
    esa estructura y esa malla en esta misma máquina.
    """
    p = Path(db_path)
    if not p.exists() or not formula or not malla:
        return 0
    kg = "x".join(str(int(x)) for x in malla)
    con = sqlite3.connect(str(p))
    try:
        fila = con.execute(
            "SELECT nk FROM calculos WHERE formula = ? AND kgrid = ? "
            "AND nk IS NOT NULL AND nk > 0 LIMIT 1", (formula, kg)).fetchone()
        return int(fila[0]) if fila else 0
    except sqlite3.Error:
        return 0
    finally:
        con.close()


# Peso de cada residuo en el ajuste: t^(-EXP_PESO).
#
# Con peso 0 (mínimos cuadrados normales) los cálculos caros dominan el
# ajuste y los baratos salen mal; con peso 1 (error relativo) pasa lo
# contrario. Se probaron los tres sobre 63 cálculos reales de 14 sistemas
# distintos, midiendo fuera de muestra: 0.5 da la menor dispersión (x1.5
# frente a x1.8) y el peor caso más pequeño (x2.6 frente a x3.6).
EXP_PESO = 0.5


def _clave_sistema(f: dict) -> tuple:
    return (f.get("natoms"), f.get("ecutwfc"), f.get("nk"),
            f.get("calculation"), f.get("nspin"))


def _prepara(filas: list) -> tuple:
    """De filas de la base a (W1, W2, T, n_scf mediana, pasos iónicos)."""
    W1, W2, T, n_scfs, pasos = [], [], [], [], {}
    for f in filas:
        d = _fila_a_descriptor(f)
        n_scf = f.get("n_scf") or N_SCF_DEFECTO
        n_ion = max(1, int(f.get("n_bfgs") or 1))
        w1, w2 = trabajo(d)
        if w1 <= 0 or not f.get("wall_s"):
            continue
        it = n_scf * n_ion
        W1.append(w1 * it); W2.append(w2 * it); T.append(float(f["wall_s"]))
        n_scfs.append(float(n_scf))
        pasos.setdefault(d["calculation"], []).append(float(n_ion))
    P = dict(PASOS_IONICOS)
    for calc, vals in pasos.items():
        P[calc] = float(statistics.median(vals))
    return (np.array(W1), np.array(W2), np.array(T),
            statistics.median(n_scfs) if n_scfs else N_SCF_DEFECTO, P)


def _ajusta(W1, W2, T) -> tuple:
    """(t0, C1, C2) por mínimos cuadrados no negativos con peso relativo."""
    if len(W1) >= 8 and W1.min() > 0 and W1.max() / W1.min() >= 5.0:
        w = np.power(np.maximum(T, 1e-9), -EXP_PESO)
        A = np.column_stack([np.ones_like(W1), W1, W2]) * w[:, None]
        sol = _nnls(A, T * w)
        if sol is not None and (sol[1] > 0 or sol[2] > 0):
            return float(sol[0]), float(sol[1]), float(sol[2])
    C = float(np.exp(np.median(np.log(T / np.maximum(W1, 1e-12)))))
    return 0.0, C, 0.0


def calibrar(db_path="olla-dft.db") -> Modelo:
    """Ajusta el coste con el historial propio y mide cuánto se equivoca.

    Se ajusta  t = t0 + C1·w1 + C2·w2  sobre los cálculos ya indexados,
    donde w1 y w2 son los dos términos de `trabajo()` por las iteraciones
    que de verdad hicieron. Los tres coeficientes significan algo: t0 es
    arrancar, C1 las FFT y C2 la ortogonalización.

    La dispersión que se devuelve NO es el residuo del ajuste. Se mide
    dejando fuera un sistema entero cada vez y prediciéndolo con los demás,
    que es la pregunta que de verdad se hace el usuario: cuánto me voy a
    equivocar con un cálculo que aún no he hecho. El residuo del ajuste da
    siempre un número más bonito y menos cierto (x1.8 frente a x1.5 sobre
    los mismos datos).
    """
    filas = historial(db_path)
    m = Modelo(fuente=str(db_path))
    if not filas:
        return m
    W1, W2, T, n_scf, P = _prepara(filas)
    if not len(W1):
        return m
    m.n = len(W1)
    m.rango = float(W1.max() / W1.min()) if W1.min() > 0 else 1.0
    m.n_scf_mediana = n_scf
    m.pasos_ionicos = P
    m.t0, m.C1, m.C2 = _ajusta(W1, W2, T)

    # --- cuánto se equivoca, medido fuera de muestra ---
    sistemas = {}
    for f in filas:
        sistemas.setdefault(_clave_sistema(f), []).append(f)
    m.n_sistemas = len(sistemas)
    razones = []
    if m.n_sistemas >= 4:
        for clave, propias in sistemas.items():
            resto = [f for k, v in sistemas.items() if k != clave for f in v]
            if len(resto) < 8:
                continue
            w1r, w2r, tr, nscf_r, Pr = _prepara(resto)
            if not len(w1r):
                continue
            t0, C1, C2 = _ajusta(w1r, w2r, tr)
            d = _fila_a_descriptor(propias[0])
            a1, a2 = trabajo(d)
            it = nscf_r * max(1.0, Pr.get(d["calculation"], 1.0))
            pred = t0 + C1 * a1 * it + C2 * a2 * it
            real = statistics.median([f["wall_s"] for f in propias])
            if pred > 0 and real > 0:
                razones.append(math.log(pred / real))
    if len(razones) >= 3:
        m.validado = True
        m.sesgo = math.exp(statistics.median(razones))
        m.dispersion = math.exp(statistics.stdev(razones))
    else:
        pred = m.t0 + m.C1 * W1 + (m.C2 or 0.0) * W2
        r = np.log(T / np.maximum(pred, 1e-9))
        m.dispersion = float(np.exp(np.std(r))) if len(r) > 1 else None
        m.sesgo = float(np.exp(np.median(r))) if len(r) else None
    return m


def _nnls(A, b):
    """Mínimos cuadrados con coeficientes no negativos.

    Un coeficiente negativo aquí no significaría nada: querría decir que
    añadir bandas hace el cálculo más rápido. Con pocos datos el ajuste
    libre los produce a menudo, así que se restringen.
    """
    try:
        from scipy.optimize import nnls
        sol, _ = nnls(A, b)
        return sol
    except Exception:                                       # noqa: BLE001
        sol, *_ = np.linalg.lstsq(A, b, rcond=None)
        return np.maximum(sol, 0.0)


# ----------------------------------------------------------------------
# Estimación
# ----------------------------------------------------------------------
def estimar(d: dict, modelo: Modelo) -> dict:
    """Segundos previstos para UN cálculo. Sin calibrar, devuelve None."""
    w1, w2 = trabajo(d)
    it = iteraciones(d, modelo)
    calc = d.get("calculation", "scf")
    fuera = {"w1": w1, "w2": w2, "iteraciones": it,
             "n_scf": modelo.n_scf_mediana if modelo.calibrado else N_SCF_DEFECTO,
             "pasos_ionicos": (modelo.pasos_ionicos.get(calc)
                               if modelo.calibrado
                               else PASOS_IONICOS.get(calc, 1.0)),
             "segundos": None, "bajo": None, "alto": None}
    if not modelo.calibrado:
        return fuera
    t = modelo.t0 + modelo.C1 * w1 * it + (modelo.C2 or 0.0) * w2 * it
    fuera["segundos"] = t
    if modelo.dispersion:
        fuera["bajo"] = t / modelo.dispersion
        fuera["alto"] = t * modelo.dispersion
    return fuera


def estimar_barrido(jobs: list, modelo: Modelo, paralelo: int = 1,
                    db_path=None) -> dict:
    """Suma el coste de todos los cálculos de un barrido."""
    total = 0.0
    detalle = []
    descs = []
    for job in jobs:
        entrada = Path(job.directory) / job.input_file
        if not entrada.exists():
            continue
        try:
            descs.append((job, descriptores_de_input(entrada)))
        except Exception:                                   # noqa: BLE001
            continue
    # Todos los puntos de un barrido comparten estructura y malla. Si UNO ya
    # corrió, su número real de puntos k vale para todos y sustituye a la
    # estimación por simetría, que es la parte más floja del modelo.
    if db_path:
        for _, d in descs:
            if d.get("nk_fuente", "").startswith(("the one", "el que")):
                continue
            n = nk_de_historial(db_path, d.get("formula"), d.get("malla"))
            if n:
                d["nk"] = n
                d["nk_fuente"] = "the one pw.x used in an identical earlier calculation"
    reales = [d["nk"] for _, d in descs if d.get("nk_fuente", "").startswith(("the one", "el que"))]
    if reales:
        real = int(statistics.median(reales))
        for _, d in descs:
            if not d.get("nk_fuente", "").startswith(("the one", "el que")):
                d["nk"] = real
                d["nk_fuente"] = "the one pw.x used at another point of the sweep"
    for job, d in descs:
        e = estimar(d, modelo)
        detalle.append((job, d, e))
        if e["segundos"]:
            total += e["segundos"]
    par = max(1, int(paralelo or 1))
    return {"detalle": detalle, "serie_s": total,
            "pared_s": total / par if total else None, "paralelo": par}


def humano(segundos) -> str:
    if segundos is None:
        return "?"
    s = float(segundos)
    if s < 90:
        return f"{s:.0f} s"
    if s < 5400:
        return f"{s / 60:.0f} min"
    if s < 172800:
        return f"{s / 3600:.1f} h"
    return f"{s / 86400:.1f} days"


def report(est: dict, modelo: Modelo) -> str:
    det = est["detalle"]
    L = ["--- Estimated cost ---"]
    if not det:
        return "Could not read the sweep inputs to estimate its cost."
    if not modelo.calibrado:
        n_at = det[0][1].get("nat", 0)
        L += [f"{len(det)} calculations, the first one with {n_at} atoms, "
              f"{det[0][1].get('nk', 1)} irreducible k-points and "
              f"{det[0][1].get('npw', 0):.0f} plane waves.",
              "",
              "Nothing to calibrate with: the calculation database is empty or does not "
              "store timings.",
              "  The model knows the SHAPE of the cost but not the speed of "
              "this machine.",
              "  Run a small sweep, index it with  olla-dft db FOLDER  and "
              "ask again:",
              "  from then on the estimates are yours, not from a generic "
              "table."]
        return "\n".join(L)

    peor = max(det, key=lambda t: t[2]["segundos"] or 0)
    fuentes = {d.get("nk_fuente", "?") for _, d, _ in det}
    det[0][1].get("calculation", "scf")
    ion = det[0][2].get("pasos_ionicos") or 1.0
    L += [f"{len(det)} calculations.  Calibrated with {modelo.n} of your earlier "
          f"calculations.",
          f"k-points: {det[0][1].get('nk', 1)} ({'; '.join(sorted(fuentes))}).",
          f"Assuming {modelo.n_scf_mediana:.0f} scf iterations"
          + (f" for each of {ion:.0f} ionic steps" if ion > 1 else "")
          + " (medians of your history).",
          ""]
    L.append(f"  {'total time in series':<26s} {humano(est['serie_s']):>10s}")
    if est["paralelo"] > 1:
        L.append(f"  {'with ' + str(est['paralelo']) + ' at a time':<26s} "
                 f"{humano(est['pared_s']):>10s}")
    L.append(f"  {'the most expensive':<26s} "
             f"{humano(peor[2]['segundos']):>10s}   ({peor[0].name})")
    if modelo.dispersion:
        lo = est["serie_s"] / modelo.dispersion / est["paralelo"]
        hi = est["serie_s"] * modelo.dispersion / est["paralelo"]
        origen = ("measured by leaving out each of your "
                  f"{modelo.n_sistemas} systems and predicting it with the rest"
                  if modelo.validado else
                  "residual of the fit, which is optimistic: there are not enough "
                  "distinct systems\n  to validate it by leaving one out")
        L += ["",
              f"Reasonable range: from {humano(lo)} to {humano(hi)}  "
              f"(spread x{modelo.dispersion:.1f}).",
              f"  {origen}.",
              "  This tool tells ten minutes from six hours, which is "
              "the decision that\n  is actually made. It is not a stopwatch."]
    if not modelo.extrapola_bien:
        L.append("\nYour history has little variety (the indexed calculations "
                 "resemble each other: the most\n  expensive does "
                 f"{modelo.rango:.0f} times more work than the cheapest). The "
                 "constant is well fitted,\n  but predicting a system "
                 "of another size may be off by a factor of two. Index\n"
                 "  calculations of different sizes and this improves by itself.")
    if any("spglib" in f for f in fuentes):
        L.append("\nThe k-points are a symmetry estimate. pw.x does not "
                 "always find the\n  same operations as spglib, and "
                 "a factor of two or three goes there. As soon as\n  you run the "
                 "first point, the estimate for the rest uses its real number.")
    if modelo.n < 8:
        L.append(f"\nThere are only {modelo.n} calculations in the database: the calibration "
                 "is weak. With fifteen or twenty\n  indexed calculations the "
                 "estimate narrows considerably.")
    return "\n".join(L)


def report_modelo(m: Modelo) -> str:
    """Qué sabe el modelo de esta máquina y cuánto se equivoca."""
    L = ["--- Cost model ---", f"Database: {m.fuente}"]
    if not m.calibrado:
        return "\n".join(L + [
            "",
            "Not calibrated: the database has no calculations with wall time.",
            "  Index folders that have already run with  olla-dft db FOLDER...  and come back.",
        ])
    L += [f"Calibrated with {m.n} calculations from {m.n_sistemas} distinct systems.",
          f"The most expensive does {m.rango:.0f} times more work than the cheapest.",
          "",
          "Coefficients (t = t0 + C1·w1 + C2·w2, per iteration):",
          f"  t0 = {m.t0:8.2f} s      startup, reading pseudos, symmetry, writing",
          f"  C1 = {m.C1:8.2e}      FFT, one per band",
          f"  C2 = {m.C2 or 0.0:8.2e}      orthogonalization, scales as bands²",
          "",
          f"Typical scf iterations: {m.n_scf_mediana:.0f}"]
    ion = {k: v for k, v in sorted(m.pasos_ionicos.items()) if v > 1}
    if ion:
        L.append("Learned ionic steps: "
                 + ", ".join(f"{k} = {v:.0f}" for k, v in ion.items()))
    L += ["", "Accuracy:"]
    if m.validado:
        L += [f"  spread  x{m.dispersion:.2f}   bias  x{m.sesgo:.2f}",
              f"  Measured by leaving out each of the {m.n_sistemas} systems "
              "and predicting it with\n  the rest. It is the real question: how much "
              "will I be off on something I have not done yet."]
        if m.sesgo and abs(math.log(m.sesgo)) > 0.15:
            direccion = "short" if m.sesgo < 1 else "long"
            L.append(f"  The model falls {direccion} on average by "
                     f"{abs(1 - m.sesgo) * 100:.0f} %. It is the normal pessimism "
                     f"of leaving\n  a system out; with more variety in the "
                     f"database it shrinks.")
    else:
        L += [f"  spread  x{m.dispersion:.2f}  (residual of the fit)",
              f"  There are only {m.n_sistemas} distinct system(s): not enough to "
              "validate by leaving one out,\n  so this number is "
              "optimistic. Index more varied calculations."]
    if not m.extrapola_bien:
        L.append("\n  Warning: the history has little variety. The scale is "
                 "well fitted, but\n  predicting a system of another size "
                 "may be off by a factor of two.")
    return "\n".join(L)
