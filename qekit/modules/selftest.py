# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Comprobación contra la física conocida, no contra uno mismo.

Las pruebas de pytest miran que el código haga lo que el código dice. El
barrido de la interfaz mira que los comandos devuelvan lo que prometen. Ni
uno ni otro detectan que una fórmula esté mal: los dos comparan Olla-DFT
consigo mismo.

Esto compara con el mundo. Cada prueba de aquí calcula una magnitud que
alguien ha medido o publicado, y la contrasta con ese valor y su fuente. Si
mañana un cambio invierte el signo del esfuerzo, o se cuela un factor 2 en
una constante, aquí sale y en las otras dos no.

Dos niveles:

  --quick   solo lo que no necesita Quantum ESPRESSO: constantes,
            integrales y límites analíticos. Segundos.
  --full    además las que corren pw.x de verdad sobre sistemas pequeños.
            Del orden de diez minutos, y hace falta un pw.x que funcione.

La comprobación con un potencial aprendido es independiente (``--mlip``):
no se ejecuta de forma implícita con ``--full`` porque MACE es opcional y no
forma parte de una instalación normal de Olla-DFT.

Las referencias llevan su procedencia. Un valor "de la literatura" sin
decir de dónde sale no es mejor que uno inventado.
"""

import shutil
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from qekit.core.errors import ErrorDeUso


@dataclass
class Prueba:
    clave: str
    titulo: str
    magnitud: str
    referencia: float
    unidad: str
    tolerancia: float           # relativa; 0.15 = 15 %
    fuente: str
    necesita_qe: bool = False
    # necesita un potencial aprendido (torch + mace): no es "rápida" aunque
    # no toque pw.x, así que se queda fuera del juego rápido
    necesita_mlip: bool = False
    coste: str = "instant"
    fn: object = None


@dataclass
class Resultado:
    prueba: Prueba
    valor: float = None
    segundos: float = 0.0
    error: str = ""
    saltada: bool = False

    @property
    def desviacion(self):
        """Desviación RELATIVA, salvo que la referencia sea cero.

        Con referencia 0 (las pruebas de invariancia: "esto tiene que dar
        exactamente cero") no hay desviación relativa que calcular, y
        dividir entre cero daba un None que luego reventaba al formatear.
        Ahí la desviación ES el valor absoluto.
        """
        if self.valor is None:
            return None
        ref = self.prueba.referencia
        if ref == 0:
            return abs(self.valor)
        return abs(self.valor - ref) / abs(ref)

    @property
    def relativa(self) -> bool:
        return self.prueba.referencia != 0

    @property
    def bien(self) -> bool:
        d = self.desviacion
        return d is not None and not self.error and d <= self.prueba.tolerancia


PRUEBAS = []


def prueba(**kw):
    def deco(fn):
        PRUEBAS.append(Prueba(fn=fn, **kw))
        return fn
    return deco


# ----------------------------------------------------------------------
# Sin Quantum ESPRESSO: constantes, integrales y límites analíticos
# ----------------------------------------------------------------------
@prueba(clave="madelung", titulo="Madelung constant, simple cubic lattice",
        magnitud="α_M", referencia=2.8372974, unidad="",
        tolerancia=1e-5,
        fuente="classic value of the Ewald sum for a point charge in "
               "a neutralizing background")
def _madelung(ctx):
    from qekit.modules import defects
    return defects.constante_madelung(np.eye(3) * 7.3)


@prueba(clave="lorenz", titulo="Lorenz number of a free-electron gas",
        magnitud="L/L₀", referencia=1.0, unidad="",
        tolerancia=0.12,
        fuente="Sommerfeld limit, L₀ = (π²/3)(k_B/e)² = 2.44e-8 W·Ω/K²")
def _lorenz(ctx):
    from qekit.modules import transport as tr
    HBAR, ME = 6.582119569e-16, 0.51099895e6 / (2.99792458e8) ** 2
    n, a = 20, 4.0
    recip = 2 * np.pi * np.linalg.inv(np.eye(3) * a).T
    ks = (np.arange(n) + 0.5) / n - 0.5
    K = np.array(np.meshgrid(ks, ks, ks, indexing="ij")).reshape(3, -1).T
    kSI = (K @ recip) * 1e10
    E = HBAR ** 2 * np.sum(kSI ** 2, axis=1) / (2 * ME)
    run = tr.TransportRun(volume=a ** 3, nelec=2.0, fermi=float(np.median(E)),
                          grid=(n, n, n))
    run.energies, run.velocities = E[:, None], (HBAR * kSI / ME)[:, None, :]
    run.weights = np.full(len(E), 2.0 / len(E))
    run = tr.compute(run, T=[300.0], mu=np.array([run.fermi]))
    return float(tr.lorenz(run, 0)[0]) / tr.L0_SOMMERFELD


@prueba(clave="npw", titulo="Plane waves of Si at 30 Ry",
        magnitud="N_PW", referencia=725.0, unidad="",
        tolerancia=0.06,
        fuente="what pw.x reports for the primitive cell of Si "
               "(V = 39.5 Å³) at 30 Ry")
def _npw(ctx):
    from qekit.modules import cost
    return cost.n_ondas_planas(39.53, 30.0)


@prueba(clave="sackur", titulo="Translational entropy of N₂ at 298 K",
        magnitud="S_trans", referencia=150.4, unidad="J/(mol·K)",
        tolerancia=0.01,
        fuente="Sackur-Tetrode at 1 bar; NIST-JANAF tables give 150.4 for the "
               "translational term")
def _sackur(ctx):
    from qekit.modules import thermochem
    s_ev = thermochem.S_traslacional(28.0134, 298.15, 100000.0)
    return s_ev * 96485.33212        # eV/(K·partícula) -> J/(mol·K)


@prueba(clave="allen_dynes", titulo="Allen-Dynes Tc for aluminium",
        magnitud="Tc", referencia=1.18, unidad="K",
        tolerancia=0.12,
        fuente="experimental Tc of Al = 1.18 K, with λ = 0.44 and "
               "ω_log = 270 K (Allen-Dynes 1975) and µ* = 0.12. NOTE: µ* is a "
               "FITTED parameter, not a computed one; with µ* = 0.10 the same "
               "formula gives 1.9 K and is not wrong for that")
def _allen(ctx):
    from qekit.modules import elph
    return elph.allen_dynes(0.44, 270.0, mustar=0.12)


@prueba(clave="allen_dynes_mu",
        titulo="Sensitivity of Allen-Dynes to µ*",
        magnitud="Tc(µ*=0.10) / Tc(µ*=0.12)", referencia=1.56, unidad="",
        tolerancia=0.05,
        fuente="the formula is exponential in µ*: raising it from 0.10 to 0.12 "
               "lowers Tc to two thirds. It is checked so that nobody cites a "
               "Tc without saying which µ* was used")
def _allen_mu(ctx):
    from qekit.modules import elph
    return (elph.allen_dynes(0.44, 270.0, mustar=0.10)
            / elph.allen_dynes(0.44, 270.0, mustar=0.12))


@prueba(clave="born2d", titulo="Sheet moduli of an isotropic sheet",
        magnitud="Y_2D", referencia=341.8, unidad="N/m",
        tolerancia=0.01,
        fuente="with C11 = 352 and C12 = 60 N/m (graphene, DFT), "
               "Y = C11 − C12²/C11")
def _born2d(ctx):
    from qekit.modules import elastic
    c11, c12 = 352.0, 60.0
    C2 = np.array([[c11, c12, 0.0], [c12, c11, 0.0],
                   [0.0, 0.0, (c11 - c12) / 2]])
    return elastic.modulos_2d(C2)["Y_x"]


@prueba(clave="gap_invariante",
        titulo="The alignment removes the arbitrary zero",
        magnitud="ΔE_v of a material with itself", referencia=0.0,
        unidad="eV", tolerancia=1e-9,
        fuente="exact identity: two calculations of the same material cannot "
               "have an offset")
def _alineamiento(ctx):
    from qekit.modules import align
    a = align.Lado(nombre="A", vbm=-5.81, cbm=-1.18, gap=4.63, referencia=0.0)
    b = align.Lado(nombre="B", vbm=131.59, cbm=136.22, gap=4.63,
                   referencia=137.40)
    return abs(align.alinear(a, b).delta_v)


@prueba(clave="ewald_escala",
        titulo="The Madelung constant does not depend on the scale",
        magnitud="|α(L=3) − α(L=30)|", referencia=0.0, unidad="",
        tolerancia=1e-6,
        fuente="exact invariance of the Ewald sum under a change of "
               "units")
def _ewald_escala(ctx):
    from qekit.modules import defects
    return abs(defects.constante_madelung(np.eye(3) * 3.0)
               - defects.constante_madelung(np.eye(3) * 30.0))


@prueba(clave="chern_qwz",
        titulo="Chern number of the Qi-Wu-Zhang insulator",
        magnitud="C (lower band, m=-1)", referencia=-1.0, unidad="",
        tolerancia=1e-10,
        fuente="Qi, Wu and Zhang, Phys. Rev. B 74, 085308 (2006): the "
               "two-band model is in a |C|=1 phase for -2<m<0; the sign "
               "here fixes the kx,ky orientation of Olla-DFT")
def _chern_qwz(ctx):
    from qekit.modules import topology

    sx = np.array([[0, 1], [1, 0]], complex)
    sy = np.array([[0, -1j], [1j, 0]], complex)
    sz = np.array([[1, 0], [0, -1]], complex)
    n = 24
    occupied = np.empty((n, n, 2, 1), complex)
    for i in range(n):
        for j in range(n):
            x, y = 2 * np.pi * i / n, 2 * np.pi * j / n
            h = (np.sin(x) * sx + np.sin(y) * sy
                 + (-1.0 + np.cos(x) + np.cos(y)) * sz)
            _energy, vectors = np.linalg.eigh(h)
            occupied[i, j] = vectors[:, :1]
    _flux, chern, _wilson, _overlap = topology.invariants_from_vectors(
        occupied)
    return chern


@prueba(clave="umklapp",
        titulo="κ_L of silicon decays as 1/T",
        magnitud="n in κ ∝ T^−n", referencia=1.0, unidad="",
        tolerancia=0.25, necesita_mlip=True, coste="~25 s",
        fuente="above the Debye temperature the phonon population "
               "grows as T and Umklapp processes are proportional "
               "to it, so κ ∝ 1/T. It is a law, not a fit, and does not "
               "depend on how good the potential is: that is why the "
               "EXPONENT is checked and not the value of κ")
def _umklapp(ctx):
    from ase.build import bulk
    from qekit.modules import kappa as kp

    si = bulk("Si", "diamond", 5.43)
    ph = kp.preparar(si, (2, 2, 2))
    s3, _ = kp.configuraciones(ph)
    F = kp.fuerzas_mlip(s3, "mace", verbose=False)
    tc, m = kp.resolver(ph, F, malla=11,
                        temperaturas=[300, 400, 500, 600, 700, 800])
    run = kp.KappaRun(formula="Si2", dim=(2, 2, 2), fuente="MACE")
    kp.recoger(run, ph, tc, m)
    return kp.exponente_temperatura(run, T_min=300.0)


@prueba(clave="her_pt", titulo="HER: platinum sits at the top of the volcano",
        magnitud="ΔG_H*", referencia=-0.09, unidad="eV",
        tolerancia=0.05,
        fuente="Nørskov et al. 2005: Pt(111) has ΔG_H* = −0.09 eV, which is "
               "why it is the best HER catalyst")
def _her_pt(ctx):
    from qekit.modules import echem
    return echem.her(-0.33).dG_H


@prueba(clave="oer_ruo2", titulo="OER: overpotential of RuO₂(110)",
        magnitud="η", referencia=0.48, unidad="V",
        tolerancia=0.10,
        fuente="Man et al. 2011 (ChemCatChem): with ΔG(OH) = 0.77, "
               "ΔG(O) = 2.16 and ΔG(OOH) = 3.87 eV, RuO₂(110) gives η ≈ 0.48 V")
def _oer_ruo2(ctx):
    from qekit.modules import echem
    e = echem.oer({"OH": 0.77, "O": 2.16, "OOH": 3.87},
                  correcciones={"OH": 0, "O": 0, "OOH": 0})
    return e.sobrepotencial


@prueba(clave="escala_oer", titulo="OOH−OH scaling relation of the OER",
        magnitud="ΔG(OOH) − ΔG(OH)", referencia=3.2, unidad="eV",
        tolerancia=0.10,
        fuente="the universal scaling relation is 3.2 ± 0.2 eV on almost "
               "every oxide surface, and from it comes the ~0.37 V limit "
               "of the OER overpotential")
def _escala(ctx):
    from qekit.modules import echem
    e = echem.oer({"OH": 0.77, "O": 2.16, "OOH": 3.87},
                  correcciones={"OH": 0, "O": 0, "OOH": 0})
    return echem.escala_ooh_oh(e)


@prueba(clave="escala_eta_min", titulo="Scaling limit of the OER "
        "overpotential",
        magnitud="η_min", referencia=0.37, unidad="V", tolerancia=0.02,
        fuente="with ΔG(OOH) − ΔG(OH) = 3.2 eV, the worst of the steps "
               "OH*→O*→OOH* is no lower than 1.6 eV, and against 1.23 V that leaves "
               "~0.37 V (Man et al. 2011)")
def _escala_eta_min(ctx):
    from qekit.modules import echem
    return echem.sobrepotencial_minimo_escala()


# ----------------------------------------------------------------------
# Con Quantum ESPRESSO
# ----------------------------------------------------------------------
def _pseudo(ctx, elemento, nombre):
    p = Path(ctx["pseudo_dir"]) / nombre
    if not p.exists():
        raise ErrorDeUso(f"{nombre} is missing in {ctx['pseudo_dir']}")
    return f"{elemento}={nombre}"


@prueba(clave="fonon_si", titulo="Optical mode of Si at Γ",
        magnitud="ω(Γ)", referencia=520.0, unidad="cm⁻¹",
        tolerancia=0.10, necesita_qe=True, coste="~20 s",
        fuente="experimental Raman of silicon: 520.7 cm⁻¹ at 300 K")
def _fonon_si(ctx):
    from ase.build import bulk
    from ase.io import write
    from qekit.modules import phonons

    d = ctx["dir"] / "fonon_si"
    d.mkdir(parents=True, exist_ok=True)
    write(d / "si.cif", bulk("Si", "diamond", 5.4073))
    run, _ = phonons.prepare(
        bulk("Si", "diamond", 5.4073), outdir=str(d), gamma_only=True,
        pseudo_dir=ctx["pseudo_dir"], ecutwfc=20, ecutrho=80, kspacing=0.4,
        insulator=True)
    ctx["correr_scf"](run.jobs, d)
    phonons.run_chain(run, pw_cmd=ctx["pw_cmd"], nproc=1, verbose=False)
    phonons.collect(run)
    f = [g[0] for g in run.gamma_freqs]
    return max(f)


@prueba(clave="wannier_si",
        titulo="Wannier centre of the Si–Si bond",
        magnitud="|r̄|", referencia=1.17563, unidad="Å",
        tolerancia=0.02, necesita_qe=True, coste="~30 s",
        fuente="the bond centre of the diamond structure is at "
               "√3·a/8 from the atom; with a = 5.43 Å that is 1.17563 Å. It is pure "
               "geometry: if the Berry phase is right, that must come out")
def _wannier_si(ctx):
    from ase.build import bulk
    from qekit.modules import wannier as wn

    d = ctx["dir"] / "wannier_si"
    si = bulk("Si", "diamond", 5.43)
    centros = [(0.125, 0.125, 0.125), (-0.375, 0.125, 0.125),
               (0.125, -0.375, 0.125), (0.125, 0.125, -0.375)]
    proy = ";".join(f"f={c[0]},{c[1]},{c[2]}:s" for c in centros)
    wn.prepare(si, outdir=str(d), malla=(4, 4, 4), proy=proy, nbnd=8,
               excluir=(5, 6, 7, 8), pseudo_dir=ctx["pseudo_dir"],
               ecutwfc=20, ecutrho=80, insulator=True)
    wn.correr(str(d), pw_cmd=ctx["pw_cmd"], nproc=1, con_bandas=False,
              verbose=False)
    run = wn.collect(str(d), atoms=si)
    if run.error_malla > wn.TOL_EXACTA:
        raise RuntimeError(
            f"the interpolation does not reproduce its own grid "
            f"({run.error_malla:.1e} eV): there is an indexing error")
    return float(np.mean([np.linalg.norm(c) for c in run.disp.centros]))


@prueba(clave="condensador",
        titulo="Charged ESM: 1/C versus distance gives 1/ε₀",
        magnitud="measured slope / (1/ε₀)", referencia=1.0, unidad="",
        tolerancia=0.06, necesita_qe=True, coste="~90 s",
        fuente="pure electrostatics: for a parallel-plate capacitor 1/C = d/ε₀. "
               "The slope does NOT depend on the material, the pseudopotential or "
               "the functional, so if the capacitance reported by Olla-DFT "
               "is a real capacitance, it must give 1/ε₀ and nothing "
               "else. It validates the formula, the area and the unit "
               "conversion at once")
def _condensador(ctx):
    from ase.build import fcc111
    from qekit.modules import esm as em

    eps0 = 8.8541878128e-12                       # F/m
    pend_ideal = 1e-10 / eps0 * 1e-2              # cm²/µF por Å
    d, invC = [], []
    for vac in (8.0, 14.0, 20.0):
        sl = fcc111("Al", size=(1, 1, 5), a=4.05, vacuum=vac / 2.0,
                    periodic=True)
        sl.center(axis=2)
        carpeta = ctx["dir"] / f"cap{int(vac)}"
        run, _c, _r = em.prepare(
            sl, outdir=str(carpeta), bc="bc3", cargas=[-0.02, 0.02],
            pseudo_dir=ctx["pseudo_dir"], ecutwfc=25, ecutrho=100,
            kspacing=0.30)
        ctx["correr"](run.jobs)
        em.collect(run, str(carpeta))
        C, _r2 = em.capacitancia(run.cargas, run.vac, run.area)
        d.append(vac / 2.0)
        invC.append(1.0 / abs(C))
    pend = float(np.polyfit(d, invC, 1)[0])
    return pend / pend_ideal


@prueba(clave="born_si", titulo="Born effective charge of silicon",
        magnitud="Z*", referencia=0.0, unidad="e",
        tolerancia=0.05, necesita_qe=True, coste="~60 s",
        fuente="in a homopolar crystal Z* is exactly ZERO, by "
               "symmetry: the two atoms are equivalent and the acoustic sum "
               "rule forces Z*₁ + Z*₂ = 0. Getting zero requires the "
               "ionic and electronic parts to cancel, and each one "
               "separately moves by 0.2 when the atom is displaced")
def _born_si(ctx):
    from ase.build import bulk
    from qekit.modules import berry as bp

    d = ctx["dir"] / "born_si"
    si = bulk("Si", "diamond", 5.43)
    run, _c, _r = bp.prepare(
        si, outdir=str(d), gdir=3, nppstr=7, kperp=(4, 4),
        desplazar=(1, [0.0, 0.0, 0.16]), nlambda=3,
        pseudo_dir=ctx["pseudo_dir"], ecutwfc=20, ecutrho=80)
    bp.correr(run, pw_cmd=ctx["pw_cmd"], nproc=1, verbose=False)
    bp.collect(run, str(d))
    comp = bp.comprobar_ionica(run)
    if comp and max(c[2] for c in comp) > 1e-4:
        raise RuntimeError(
            "the ionic phase does not match Σ Z_a·f_a: the geometry or the "
            "valences are not what the module believes")
    return bp.analizar(run)["zeff"]


@prueba(clave="gamma_al", titulo="Surface energy of Al(111)",
        magnitud="γ", referencia=1.10, unidad="J/m²",
        tolerancia=0.25, necesita_qe=True, coste="~60 s",
        fuente="full-potential LDA (Vitos et al. 1998) gives 1.20 J/m²; "
               "the polycrystalline experiment, 1.14")
def _gamma_al(ctx):
    from ase.build import bulk
    from qekit.modules import surfen

    d = ctx["dir"] / "gamma_al"
    run, _ = surfen.prepare(
        bulk("Al", "fcc", 4.05), miller=(1, 1, 1), capas=(4, 5, 6, 7),
        vacuum=16.0, outdir=str(d), pseudo_dir=ctx["pseudo_dir"],
        ecutwfc=24, ecutrho=96, kspacing=0.25)
    ctx["correr"](run.jobs)
    surfen.collect(run)
    if run.gamma_ajuste is None:
        raise RuntimeError("the fit did not work out")
    return run.gamma_ajuste * surfen.EV_A2_A_J_M2


@prueba(clave="bulk_si", titulo="Bulk modulus of Si by strain",
        magnitud="B", referencia=95.0, unidad="GPa",
        tolerancia=0.15, necesita_qe=True, coste="~50 s",
        fuente="LDA gives 93-97 GPa (Nielsen & Martin 1985); the experiment, 98")
def _bulk_si(ctx):
    from ase.build import bulk
    from qekit.modules import strain

    d = ctx["dir"] / "bulk_si"
    run, _ = strain.prepare(
        bulk("Si", "diamond", 5.4073), modo="hidrostatica", rangos="-2:2:5",
        outdir=str(d), pseudo_dir=ctx["pseudo_dir"], ecutwfc=32, ecutrho=128,
        kspacing=0.25, relax_ions=False, insulator=True)
    ctx["correr"](run.jobs)
    strain.collect(run)
    idx = [i for i in run.ok if run.energies[i] is not None]
    if len(idx) < 3:
        raise RuntimeError("not enough points came out")
    x = np.array([run.strains[i] for i in idx])
    y = np.array([run.energies[i] for i in idx])
    a2 = np.polyfit(x, y, 2)[0]
    return 2 * a2 / (9 * run.volume0) * 160.21766


@prueba(clave="sitio_h_al", titulo="H on Al(111): the hollow beats the top",
        magnitud="E_ads(top) − E_ads(hollow)", referencia=5.6, unidad="eV",
        tolerancia=0.60, necesita_qe=True, coste="~60 s",
        fuente="hydrogen chemisorbs in the hollow of fcc(111), not on top "
               "of an atom; the order hollow < bridge < top is textbook")
def _sitio_h(ctx):
    from ase.build import fcc111
    from qekit.modules import adsorb

    d = ctx["dir"] / "ads_h"
    sl = fcc111("Al", size=(2, 2, 3), vacuum=8.0)
    sl.pbc = (True, True, True)
    run, _ = adsorb.prepare(
        sl, "H", outdir=str(d), altura=1.0, relax_ions=False,
        pseudo_dir=ctx["pseudo_dir"], ecutwfc=20, ecutrho=80, kspacing=0.30)
    ctx["correr"](run.jobs)
    adsorb.collect(run)
    e = run.energias_ads
    por_tipo = {}
    for i, s in enumerate(run.sitios):
        if e[i] is not None:
            por_tipo.setdefault(s.tipo, []).append(e[i])
    if "top" not in por_tipo or "hollow" not in por_tipo:
        raise RuntimeError("sites are missing")
    return min(por_tipo["top"]) - min(por_tipo["hollow"])


# ----------------------------------------------------------------------
# Ejecución
# ----------------------------------------------------------------------
def ejecutar(claves=None, con_qe: bool = False, con_mlip: bool = False,
             pseudo_dir: str = None, pw_cmd: str = None, nproc: int = None,
             paralelo: int = 1, carpeta: str = None,
             verbose: bool = True) -> list:
    from qekit.core import runner as run_mod

    seleccion = [
        p for p in PRUEBAS
        if (claves is None or p.clave in claves)
        and (not p.necesita_qe or con_qe)
        and (not p.necesita_mlip or con_mlip)
    ]
    if not seleccion:
        raise ErrorDeUso(
            "no test matches. The available ones: "
            + ", ".join(p.clave for p in PRUEBAS))

    tmp = Path(carpeta) if carpeta else Path(tempfile.mkdtemp(prefix="qekit_st_"))
    tmp.mkdir(parents=True, exist_ok=True)

    def _correr(jobs):
        res = run_mod.run_all(jobs, pw_cmd=pw_cmd, nproc=nproc,
                              paralelo=paralelo, verbose=False)
        malos = [r for r in res if not r.ok]
        if malos:
            raise RuntimeError(f"{len(malos)} of {len(jobs)} calculations "
                               f"failed: {malos[0].error}")
        return res

    def _correr_scf(jobs, d):
        return _correr(jobs)

    ctx = {"dir": tmp, "pseudo_dir": pseudo_dir, "pw_cmd": pw_cmd,
           "correr": _correr, "correr_scf": _correr_scf}

    fuera = []
    for p in seleccion:
        if verbose:
            print(f"  {p.clave:16s} ", end="", flush=True)
        r = Resultado(prueba=p)
        t0 = time.time()
        try:
            r.valor = float(p.fn(ctx))
        except Exception as exc:                            # noqa: BLE001
            r.error = f"{type(exc).__name__}: {exc}"
        r.segundos = time.time() - t0
        fuera.append(r)
        if verbose:
            if r.error:
                print(f"ERROR  ({r.error[:60]})")
            else:
                marca = "ok " if r.bien else "BAD"
                detalle = (f"({r.desviacion * 100:5.2f} % of {p.referencia:g})"
                           if r.relativa
                           else f"(must be 0; gives {r.desviacion:.1e})")
                print(f"{marca}  {r.valor:12.5g} {p.unidad:8s} {detalle}")
    if carpeta is None:
        shutil.rmtree(tmp, ignore_errors=True)
    return fuera


def report(resultados: list) -> str:
    bien = [r for r in resultados if r.bien]
    mal = [r for r in resultados if not r.bien and not r.error]
    err = [r for r in resultados if r.error]
    L = ["--- Validation against known physics ---",
         f"{len(resultados)} tests: {len(bien)} passed, {len(mal)} out of "
         f"tolerance, {len(err)} with error",
         ""]
    for r in resultados:
        p = r.prueba
        estado = "ERROR" if r.error else ("ok" if r.bien else "BAD")
        L.append(f"[{estado:^5s}] {p.titulo}")
        if r.error:
            L.append(f"          {r.error}")
        else:
            if r.relativa:
                L.append(f"          {p.magnitud} = {r.valor:.5g} {p.unidad}"
                         f"   reference {p.referencia:g} {p.unidad}"
                         f"   ({r.desviacion * 100:.2f} %, tolerance "
                         f"{p.tolerancia * 100:.0f} %)")
            else:
                L.append(f"          {p.magnitud} = {r.valor:.3e} {p.unidad}"
                         f"   must be 0 (tolerance "
                         f"{p.tolerancia:g})")
        L.append(f"          source: {p.fuente}")
        L.append("")
    if mal:
        L += ["Those that come out BAD are not always a code failure: a "
              "tight\n  tolerance, a different pseudopotential or a "
              "low cutoff also move them.\n  What they do mean is "
              "that this number has changed and one has to look at why.", ""]
    if err:
        L += ["Those that give ERROR did not manage to produce a number: a "
              "pseudopotential is missing,\n  pw.x is not there, or the calculation did not converge.", ""]
    return "\n".join(L)
