# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Conductividad térmica de red: el fonón que se dispersa contra otro fonón.

Los fonones armónicos no conducen calor de forma finita: si el cristal fuera
exactamente armónico, un paquete de fonones viajaría para siempre y κ sería
infinita. La conductividad sale del término CÚBICO de la energía, el que
permite que un fonón se parta en dos o que dos se fundan en uno. Por eso
hace falta la tercera derivada de la energía respecto de los
desplazamientos, la fc3, y por eso este cálculo es tan caro: la fc2 necesita
una derivada por par de átomos y la fc3, una por TRÍO.

En la práctica:

    Φ_ijk = ∂³E / ∂u_i ∂u_j ∂u_k

se obtiene por diferencias finitas, desplazando dos átomos a la vez en una
supercelda. El número de configuraciones crece deprisa con el tamaño de la
supercelda, y ahí está todo el coste. Con las fuerzas ya calculadas, la
ecuación de Boltzmann de fonones en aproximación de tiempo de relajación da

    κ_L = (1/NV) Σ_λ C_λ v_λ ⊗ v_λ τ_λ

que es lo que resuelve phono3py y lo que se reporta aquí.

Este módulo permite dos fuentes de fuerzas, y la diferencia importa:

  - **Quantum ESPRESSO.** Es el cálculo de verdad. En silicio, una supercelda
    2×2×2 de la primitiva son 57 configuraciones de 16 átomos.
  - **Un potencial aprendido (MACE y compañía).** Las mismas 57
    configuraciones en 8 segundos en vez de 40 minutos. Sirve para explorar
    y para decidir el tamaño de la supercelda antes de gastar el cálculo
    caro, pero **el valor absoluto puede estar muy lejos**: con MACE-MP
    pequeño el silicio da 51 W/mK a 300 K frente a los ~140 medidos. La
    dependencia con T sale bien (el 1/T de Umklapp), el número no. Se avisa
    en el informe, siempre.

Lo que hay que converger, y son tres cosas a la vez:

  1. El tamaño de la supercelda de la fc3. Es lo que fija hasta qué distancia
     se tienen en cuenta las interacciones anarmónicas.
  2. La malla de q sobre la que se resuelve la ecuación de Boltzmann.
  3. Y, si se compara con un experimento, los isótopos: el silicio natural
     conduce alrededor de un 10 % menos que el isotópicamente puro, y esa
     diferencia es mayor que muchas de las que se discuten en un artículo.

La aproximación de tiempo de relajación (RTA) subestima κ frente a la
solución exacta de la ecuación de Boltzmann, típicamente un 10-15 % en
silicio y mucho más en materiales con procesos normales dominantes, como el
grafeno o el diamante. No es un error: es la aproximación, y se dice.
"""

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from qekit.core.errors import FaltanDatos
from qekit.core import style as qstyle

# recorridos libres medios (Å) sobre los que se acumula κ
RECORRIDOS = np.logspace(0, 7, 141)
# por encima de esto se avisa antes de escribir cientos de inputs
MUCHAS_CONFIGURACIONES = 150


def _phonopy_atoms(atoms):
    from phonopy.structure.atoms import PhonopyAtoms
    return PhonopyAtoms(symbols=atoms.get_chemical_symbols(),
                        cell=np.array(atoms.cell.array, float),
                        scaled_positions=atoms.get_scaled_positions())


def _ase(sc):
    from ase import Atoms
    return Atoms(symbols=list(sc.symbols), positions=np.array(sc.positions),
                 cell=np.array(sc.cell), pbc=True)


def _importar():
    try:
        from phono3py import Phono3py
    except ImportError as exc:
        raise FaltanDatos(
            "phono3py is required, which is what solves the phonon "
            "Boltzmann equation:\n  pip install phono3py\n"
            "It is installed with pip and does not require recompiling Quantum "
            "ESPRESSO.") from exc
    return Phono3py


@dataclass
class KappaRun:
    """Un cálculo de conductividad térmica de red."""
    formula: str = ""
    dim: tuple = (2, 2, 2)
    dim_fc2: tuple = None
    distancia: float = 0.03
    n_config: int = 0
    n_atomos: int = 0
    fuente: str = ""            # 'Quantum ESPRESSO' o el nombre del potencial
    malla: tuple = (11, 11, 11)
    temperaturas: np.ndarray = None
    kappa: np.ndarray = None    # (nT, 6) en notación de Voigt, W/mK
    isotopos: bool = False
    frontera: float = None      # tamaño de grano en µm, si se puso
    gamma: np.ndarray = None
    frecuencias: np.ndarray = None
    pesos: np.ndarray = None
    velocidades: np.ndarray = None
    cv: np.ndarray = None
    mode_kappa: np.ndarray = None   # (nT, gp, banda, 6) κ por modo, de phono3py
    gamma_iso: np.ndarray = None    # (gp, banda) dispersión por isótopos
    i300: int = None
    avisos: list = field(default_factory=list)
    directorio: str = ""

    @property
    def kappa_media(self):
        """(κ_xx + κ_yy + κ_zz)/3, que es lo que se compara con un policristal."""
        if self.kappa is None:
            return None
        return self.kappa[:, :3].mean(axis=1)


def preparar(atoms, dim=(2, 2, 2), dim_fc2=None, distancia=0.03,
             simetria=1e-5):
    """Genera las configuraciones desplazadas de la fc3 (y de la fc2).

    `dim_fc2` permite una supercelda MAYOR solo para la parte armónica, que
    es mucho más barata y a la vez la que más necesita alcance: las
    constantes de fuerza armónicas de un semiconductor llegan lejos, las
    anarmónicas no tanto. Es la forma estándar de no pagar la fc3 en una
    supercelda enorme.
    """
    Phono3py = _importar()
    ph = Phono3py(_phonopy_atoms(atoms),
                  supercell_matrix=list(map(int, dim)),
                  phonon_supercell_matrix=(list(map(int, dim_fc2))
                                           if dim_fc2 else None),
                  primitive_matrix="auto", symprec=simetria, log_level=0)
    ph.generate_displacements(distance=float(distancia))
    return ph


def configuraciones(ph):
    """Las superceldas desplazadas, como objetos de ASE."""
    sc3 = [_ase(s) for s in ph.supercells_with_displacements if s is not None]
    sc2 = []
    if ph.phonon_supercell_matrix is not None:
        sc2 = [_ase(s) for s in ph.phonon_supercells_with_displacements
               if s is not None]
    return sc3, sc2


def fuerzas_mlip(configs, modelo="mace", verbose=True):
    """Fuerzas con un potencial aprendido. Segundos en vez de horas."""
    from qekit.modules import mlip
    calc = mlip.calculator(modelo)
    fuera = []
    for i, a in enumerate(configs):
        a = a.copy()
        a.calc = calc
        fuera.append(a.get_forces())
        if verbose and (i % 25 == 0 or i == len(configs) - 1):
            print(f"    {i + 1}/{len(configs)}", end="\r", flush=True)
    if verbose:
        print(" " * 30, end="\r")
    return np.array(fuera)


def escribir_inputs(configs, destino, common, prefijo="d", kspacing=None,
                    ecutwfc=None, ecutrho=None, conv_thr=1e-10):
    """Un scf con cálculo de fuerzas por configuración, cada uno en su carpeta."""
    from qekit.modules import inputgen, sweep

    destino = Path(destino); destino.mkdir(parents=True, exist_ok=True)
    carpetas = []
    for i, a in enumerate(configs):
        d = destino / f"{prefijo}{i:04d}"
        d.mkdir(parents=True, exist_ok=True)
        grid = sweep.default_grid(a, kspacing)
        txt = inputgen.build_pw_input(
            atoms=a, pseudos=common["pseudos"], calculation="scf",
            prefix=common["prefix"], pseudo_dir=common["pseudo_dir"],
            ecutwfc=ecutwfc or common["ecutwfc"],
            ecutrho=ecutrho or common["ecutrho"],
            kcard=f"K_POINTS automatic\n  {grid[0]} {grid[1]} {grid[2]} "
                  "0 0 0\n",
            insulator=common["insulator"], degauss=common["degauss"],
            smearing=common["smearing"], conv_thr=conv_thr)
        sweep.write_input(d / "pw.in", txt)
        carpetas.append(d)
    return carpetas


def leer_fuerzas(carpetas, natomos):
    """Fuerzas de cada scf, en eV/Å y en el orden de los átomos del input."""
    from qekit.core import qeout
    fuera, faltan = [], []
    for d in carpetas:
        d = Path(d)
        try:
            res = qeout.read_xml(str(d))
            F = np.array(res.forces, float)
        except Exception:                                   # noqa: BLE001
            faltan.append(d.name)
            continue
        if F is None or F.shape != (natomos, 3):
            faltan.append(d.name)
            continue
        fuera.append(F)
    if faltan:
        raise FaltanDatos(
            f"the forces of {len(faltan)} configurations are missing "
            f"({', '.join(faltan[:5])}{'...' if len(faltan) > 5 else ''}).\n"
            f"Without ALL of them the fc3 cannot be built: each one contributes a "
            f"different derivative.")
    return np.array(fuera)


def resolver(ph, fuerzas, fuerzas_fc2=None, malla=11, temperaturas=None,
             isotopos=False, frontera_um=None, simetrizar=True):
    """Construye fc2 y fc3 y resuelve la ecuación de Boltzmann en RTA.

    `frontera_um` añade dispersión por el tamaño de grano, que es lo que
    convierte este cálculo en algo comparable con una película delgada o un
    nanohilo. `isotopos` añade la dispersión por masa con las abundancias
    naturales: en silicio son ~10 %, más de lo que suele discutirse en un
    artículo.
    """
    temperaturas = (np.asarray(temperaturas, float) if temperaturas is not None
                    else np.arange(100.0, 801.0, 100.0))
    ph.forces = np.asarray(fuerzas, float)
    if fuerzas_fc2 is not None and len(np.asarray(fuerzas_fc2)):
        ph.phonon_forces = np.asarray(fuerzas_fc2, float)
    ph.produce_fc3()
    ph.produce_fc2()
    if simetrizar:
        ph.symmetrize_fc3()
        ph.symmetrize_fc2()
    m = (int(malla),) * 3 if np.isscalar(malla) else tuple(int(x) for x in malla)
    ph.mesh_numbers = list(m)
    ph.init_phph_interaction()
    ph.run_thermal_conductivity(
        temperatures=temperaturas,
        is_isotope=bool(isotopos),
        # boundary_mfp de phono3py está EN MICRÓMETROS: su propio código
        # calcula Γ_b = |v|·1e6·Angstrom/(4π·boundary_mfp) y lo documenta así
        # (scattering_solvers.py). Pasarlo en Å multiplicando por 1e4 hacía la
        # dispersión por fronteras 10⁴ veces más débil de lo pedido: --grain
        # quedaba inerte y κ salía la del cristal infinito, mientras el informe
        # decía que se había aplicado. El 1e6 de «sin fronteras» (= 1 m) ya
        # estaba en µm, que es lo que delataba la incoherencia.
        boundary_mfp=(float(frontera_um) if frontera_um else 1e6),
        write_kappa=False)
    return ph.thermal_conductivity, m


def recoger(run, ph, tc, malla):
    """Pasa lo que devuelve phono3py a los campos del KappaRun."""
    run.malla = tuple(malla)
    run.temperaturas = np.asarray(tc.temperatures, float)
    run.kappa = np.asarray(tc.kappa[0], float)          # (nT, 6) Voigt
    run.frecuencias = np.asarray(tc.frequencies, float)
    run.pesos = np.asarray(tc.grid_weights, float)
    # Estos tres son los que alimentan la curva de recorrido libre medio. Si
    # phono3py renombra uno, antes desaparecían EN SILENCIO la sección del
    # informe, KAPPA_recorrido.dat y la segunda figura, y el usuario solo veía
    # que faltaba salida. Ahora se acota a lo que puede pasar de verdad
    # (un atributo que ya no está) y queda dicho en los avisos.
    try:
        run.gamma = np.asarray(tc.gamma[0], float)
        run.velocidades = np.asarray(tc.group_velocities, float)
        run.cv = np.asarray(tc.mode_heat_capacities, float)
        # κ POR MODO, tal como la reparte phono3py. Es lo que hace que la
        # acumulada sume por construcción la κ que se informa arriba, en vez
        # de una reconstrucción a mano que se desviaba de ella.
        mk = getattr(tc, "mode_kappa", None)
        run.mode_kappa = None if mk is None else np.asarray(mk[0], float)
        iso = getattr(tc, "gamma_isotope", None)
        run.gamma_iso = None if iso is None else np.asarray(iso[0], float)
    except (AttributeError, TypeError, IndexError) as exc:
        run.gamma = run.velocidades = run.cv = None
        run.mode_kappa = run.gamma_iso = None
        run.avisos.append(
            "this phono3py does not expose the per-mode data (gamma, group "
            f"velocities, heat capacities): {exc}.\n  κ(T) is unaffected, but "
            "there is no mean free path section, no KAPPA_recorrido.dat and no "
            "second figure.")
    d = np.abs(run.temperaturas - 300.0)
    run.i300 = int(np.argmin(d)) if len(d) else None
    return run


# 1e-10 m: el «Angstrom» que phono3py mete en su término de fronteras.
_ANGSTROM = 1e-10


def gamma_total(run, iT):
    """Γ efectiva del modo: fonón-fonón + isótopos + fronteras.

    Es la suma que phono3py usa para κ (`grid_point_data.py`), y la que hay
    que usar para Λ. Con solo la parte fonón-fonón, `--isotopes` y `--grain`
    cambiaban la κ de la tabla pero NO el recorrido libre medio de debajo, que
    seguía siendo el del cristal puro e infinito sin decirlo.
    """
    g = np.array(run.gamma[iT], dtype=float)
    if run.gamma_iso is not None:
        g = g + run.gamma_iso
    if run.frontera and run.velocidades is not None:
        # Γ_b = |v|·1e6·Angstrom/(4π·L), con L en µm; es literalmente la
        # fórmula de phono3py (scattering_solvers.py).
        vmod = np.sqrt((run.velocidades ** 2).sum(axis=-1))
        g = g + vmod * 1e6 * _ANGSTROM / (4.0 * np.pi * float(run.frontera))
    return g


def acumulada(run, iT=None):
    """κ acumulada frente al recorrido libre medio.

    Es la curva que dice si nanoestructurar sirve: si el 80 % de κ lo llevan
    fonones con Λ > 100 nm, un grano de 50 nm corta ese 80 %. Si lo llevan
    fonones de 5 nm, no hay nada que hacer con el tamaño de grano.

    El peso de cada modo es su `mode_kappa`, la descomposición que phono3py ya
    tiene hecha, así que la curva suma EXACTAMENTE la κ que se informa arriba.
    Reconstruirla a mano como C·v²·τ/3 daba una curva que no cerraba con ella
    en cuanto entraba cualquier otro canal de dispersión. Si la librería no la
    expone (versiones viejas), se usa esa reconstrucción como respaldo.
    """
    if run.gamma is None or run.velocidades is None or run.cv is None:
        return None, None
    iT = run.i300 if iT is None else int(iT)
    if iT is None or run.pesos is None:
        return None, None
    vmod = np.sqrt((run.velocidades ** 2).sum(axis=-1))
    # Los modos acústicos en Γ tienen Γ = 0 exactamente: τ es infinito y el
    # producto τ·v² es 0·∞ = NaN. No es un fallo, es que esos modos no
    # transportan calor en un cristal infinito; se descartan más abajo. El
    # errestate tiene que cubrir TAMBIÉN los productos, no solo la división.
    with np.errstate(divide="ignore", invalid="ignore"):
        # τ = 1/(2·2π·Γ), el convenio del propio phono3py (su get_mfp calcula
        # |v|/(2·2π·Γ) y su factor a W/mK lleva el mismo 1/2π «from definition
        # of lifetime»). Un 2 pasa la Γ de HWHM a anchura total y el 2π pasa
        # de frecuencia cíclica a angular, porque la Γ que devuelve viene en
        # THz ordinarios. Sin el 2π, Λ salía 6.28 veces más largo.
        L = vmod / (2.0 * 2.0 * np.pi * gamma_total(run, iT))     # Å
        if run.mode_kappa is not None:
            # traza/3 de la κ del modo en notación de Voigt; ya lleva el peso
            # del punto q, así que NO se vuelve a multiplicar por run.pesos
            contrib = np.asarray(run.mode_kappa[iT])[..., :3].mean(axis=-1)
            w = np.ones_like(contrib)
        else:
            v2 = (run.velocidades ** 2).sum(axis=-1) / 3.0
            tau = 1.0 / (2.0 * 2.0 * np.pi * gamma_total(run, iT))
            contrib = run.cv[iT] * v2 * tau
            w = run.pesos[:, None] * np.ones_like(contrib)
    ok = np.isfinite(L) & np.isfinite(contrib) & (contrib > 0)
    Lf, cf, wf = L[ok], contrib[ok], w[ok]
    if Lf.size == 0:                    # ningún modo aporta (todo Γ o NaN)
        return None, None
    orden = np.argsort(Lf)
    acum = np.cumsum(cf[orden] * wf[orden])
    if acum[-1] <= 0:
        return None, None
    return Lf[orden], acum / acum[-1]


def recorrido_representativo(run, fraccion=0.5, iT=None):
    """El Λ por debajo del cual se acumula `fraccion` de κ."""
    L, a = acumulada(run, iT)
    if L is None:
        return None
    return float(np.interp(float(fraccion), a, L))


def exponente_temperatura(run, T_min=200.0):
    """El n de κ ∝ T^(−n). Los procesos Umklapp puros dan n = 1."""
    if run.kappa is None or run.temperaturas is None:
        return None
    k = run.kappa_media
    m = (run.temperaturas >= T_min) & (k > 0)
    if m.sum() < 3:
        return None
    return float(-np.polyfit(np.log(run.temperaturas[m]), np.log(k[m]), 1)[0])


def report(run) -> str:
    # Los avisos se ARMAN aquí, en una lista local. Antes se hacía `append`
    # sobre run.avisos, y como el CLI llama a report() y después a export()
    # —que vuelve a llamarlo—, cada WARNING acababa duplicado en KAPPA.txt, y
    # triplicado si algo más lo llamaba una tercera vez.
    avisos = list(run.avisos)
    L = [f"--- Lattice thermal conductivity: {run.formula} ---",
         f"Forces: {run.fuente}",
         f"fc3 supercell: {run.dim[0]}×{run.dim[1]}×{run.dim[2]}"
         + (f"   |   fc2: {run.dim_fc2[0]}×{run.dim_fc2[1]}×{run.dim_fc2[2]}"
            if run.dim_fc2 else "")
         + f"   |   {run.n_config} configurations of {run.n_atomos} atoms",
         f"q-grid: {run.malla[0]}×{run.malla[1]}×{run.malla[2]}"
         + ("   |   natural isotopes" if run.isotopos else "")
         + (f"   |   grains of {run.frontera:g} µm" if run.frontera else ""),
         ""]
    if run.kappa is None:
        return "\n".join(L + ["No κ yet."])
    iso = np.allclose(run.kappa[:, 0], run.kappa[:, 1], rtol=0.02) and \
        np.allclose(run.kappa[:, 0], run.kappa[:, 2], rtol=0.02)
    L += ["   T (K)      κ_xx      κ_yy      κ_zz      mean  (W/m·K)"]
    for i, T in enumerate(run.temperaturas):
        k = run.kappa[i]
        marca = "  ←" if i == run.i300 else ""
        L.append(f"  {T:6.0f}  {k[0]:9.2f} {k[1]:9.2f} {k[2]:9.2f} "
                 f"{run.kappa_media[i]:9.2f}{marca}")
    if iso:
        L.append("  (the tensor is isotropic, as befits cubic "
                 "symmetry)")
    n = exponente_temperatura(run)
    if n is not None:
        L += ["", f"Temperature dependence: κ ∝ T^−{n:.2f}"]
        if abs(n - 1.0) < 0.25:
            L.append("  This is the T⁻¹ of Umklapp processes: above the "
                     "Debye temperature,")
            L.append("  the more phonons there are, the more they scatter each other.")
        else:
            L.append("  It departs from the T⁻¹ of pure Umklapp; this usually means "
                     "there is another dominant")
            L.append("  channel (boundaries, isotopes) or that the q-grid "
                     "is not converged.")
    L50 = recorrido_representativo(run, 0.5)
    L90 = recorrido_representativo(run, 0.9)
    if L50:
        L += ["", "Mean free path (at the T closest to 300 K):",
              f"  half of κ is carried by phonons with Λ < {L50 / 10:.0f} nm",
              f"  90 %,                                   Λ < {L90 / 10:.0f} nm",
              "  This is what tells whether nanostructuring helps: a grain "
              "smaller than that Λ cuts",
              "  that part of κ; a larger one does nothing."]
        # Λ90 vive en la cola de fonones de recorrido largo, que es justo lo
        # que una malla q floja no muestrea. Medido en Si con Stillinger-Weber
        # (supercelda 3×3×3): de 11³ a 31³, Λ50 se mueve de 0.68 a 0.82 µm
        # (+21 %) pero Λ90 pasa de 4.4 a 21.9 µm, cinco veces. Λ50 se puede
        # citar con una malla normal; Λ90 no, y antes no se decía.
        if run.malla and min(run.malla) < 25:
            avisos.append(
                f"the {run.malla[0]}×{run.malla[1]}×{run.malla[2]} q-grid is "
                "enough for Λ50 but NOT for Λ90.\n  The 90 % figure lives in "
                "the long-mean-free-path tail, which a coarse grid does not "
                "sample:\n  in silicon it grows by a factor of five between "
                "11³ and 31³ while Λ50 moves 20 %.\n  Converge the grid "
                "before quoting Λ90 or sizing a grain from it.")
    L += ["", "What is NOT included, and is worth keeping in mind:",
          "  · This is RTA, not the exact solution of the Boltzmann equation. "
          "RTA underestimates κ",
          "    (≈10-15 % in silicon, much more in graphene or diamond).",
          "  · Only three-phonon scattering. At high temperature "
          "four-phonon processes",
          "    lower κ, and in very anharmonic materials they are not a "
          "detail."]
    if not run.isotopos:
        L.append("  · No isotope scattering. Natural silicon conducts "
                 "~10 % less than the")
        L.append("    isotopically pure one: if you compare with an experiment, "
                 "use --isotopes.")
    if run.fuente and "ESPRESSO" not in run.fuente.upper():
        avisos.append(
            f"The forces come from {run.fuente}, not from DFT. The shape of "
            f"κ(T) usually comes out right,\n  but the absolute value may be "
            f"far off: with small MACE-MP silicon gives\n  ~51 W/mK at 300 K "
            f"where the experiment is ~140. Use it to choose the supercell "
            f"and\n  the grid, and repeat with Quantum ESPRESSO before "
            f"publishing anything.")
    if int(np.prod(run.dim)) <= 8:
        avisos.append(
            f"The fc3 supercell is {run.dim[0]}×{run.dim[1]}×{run.dim[2]}"
            f", which is small. κ has to converge\n  in the supercell "
            f"size AND in the q-grid at the same time: raise one, then the "
            f"other,\n  and do not trust it until neither of the two moves the "
            f"result.")
    for a in avisos:
        L += ["", f"WARNING: {a}"]
    return "\n".join(L)


def export(run, outdir="kappa") -> list:
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    escritos = []
    if run.kappa is not None:
        f = out / "KAPPA.dat"
        np.savetxt(f, np.column_stack([run.temperaturas, run.kappa,
                                       run.kappa_media]), fmt="%12.5f",
                   header="T(K)  kxx  kyy  kzz  kyz  kxz  kxy  mean "
                          "(W/m/K)")
        escritos.append(str(f))
    L, a = acumulada(run)
    if L is not None:
        f = out / "KAPPA_recorrido.dat"
        Lg = RECORRIDOS
        ac = np.interp(Lg, L, a)
        np.savetxt(f, np.column_stack([Lg / 10.0, ac]), fmt="%14.6e",
                   header="Lambda(nm)   cumulative fraction of kappa")
        escritos.append(str(f))
    f = out / "KAPPA.txt"
    f.write_text(report(run) + "\n", encoding="utf-8")
    escritos.append(str(f))
    return escritos


def plot(run, outfile="kappa", formats="pdf,png", theme=None, size=None,
         family=None, background=None, palette=None, usetex=None,
         width="single", journal="generic", aspect=0.72, mono=False,
         dpi=None) -> list:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:                              # pragma: no cover
        raise RuntimeError("matplotlib is not installed.") from exc
    if run.kappa is None:
        raise FaltanDatos("there is no κ to plot.")
    st = qstyle.apply(theme, size=size, family=family, background=background,
                      palette=palette, usetex=usetex, mono=mono)
    fig, ax = qstyle.new_figure(width, journal, aspect)
    cols = qstyle.palette(3, mono=mono)
    ax.loglog(run.temperaturas, run.kappa_media, marker="o", ms=4,
              lw=st["line"], color=cols[0], label=r"$\kappa_L$ computed")
    n = exponente_temperatura(run)
    if n is not None:
        i = run.i300 if run.i300 is not None else 0
        T0, k0 = run.temperaturas[i], run.kappa_media[i]
        ax.loglog(run.temperaturas, k0 * (run.temperaturas / T0) ** (-1.0),
                  lw=st["line"], color=cols[1], dashes=[4.0, 2.0],
                  label=r"$T^{-1}$ (Umklapp)")
    ax.set_xlabel("T (K)")
    ax.set_ylabel(r"$\kappa_L$ (W m$^{-1}$K$^{-1}$)")
    ax.legend(frameon=False, fontsize=st["legend"])
    escritos = qstyle.save(fig, outfile, formats, dpi=dpi, modulo="kappa")
    plt.close(fig)

    L, a = acumulada(run)
    if L is not None:
        fig2, ax2 = qstyle.new_figure(width, journal, aspect)
        ax2.semilogx(L / 10.0, a * 100.0, lw=st["line"], color=cols[0])
        for frac, et in ((0.5, "50 %"), (0.9, "90 %")):
            x = recorrido_representativo(run, frac) / 10.0
            ax2.axvline(x, color=cols[1], lw=0.7, dashes=[3.0, 2.0])
            ax2.annotate(f"{et}: {x:.0f} nm", (x, frac * 100),
                         textcoords="offset points", xytext=(5, -10),
                         fontsize=st["legend"], color=cols[1])
        ax2.set_xlabel(r"mean free path $\Lambda$ (nm)")
        ax2.set_ylabel(r"cumulative % of $\kappa_L$")
        ax2.set_ylim(0, 102)
        escritos += qstyle.save(fig2, str(outfile) + "_recorrido", formats,
                                dpi=dpi, modulo="kappa")
        plt.close(fig2)
    return escritos
