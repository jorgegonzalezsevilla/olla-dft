"""Cada fórmula contra una implementación de referencia o una identidad exacta.

La auditoría de la 1.6.0 dejó 52 afirmaciones sin adjudicar porque la fuente
primaria estaba tras un muro de pago o el proxy la bloqueaba. Lo que sí cierra
ese hueco — y es lo que destapó los errores de `kappa`, del signo de la HER y
de la transformada de dE/dk — es comprobar el número contra otra cosa que ya
sabe la respuesta: ASE, o una identidad física que no depende de que nadie
publique nada.

Estas pruebas son lentas de escribir y rápidas de correr. Fijan el CONVENIO,
que es lo que se rompe sin que nadie lo note: el 1/6 del MSD, el signo del
Seebeck, el 2π de una frecuencia, qué n es la n de Slack.
"""
import numpy as np
import pytest

from qekit.core.compat import trapezoid


# ----------------------------------------------------------------------
# Termoquímica contra ase.thermochemistry
# ----------------------------------------------------------------------
# ASE trae su propia edición de CODATA (la 3.29 usa la de 2014, el paquete
# usa la de 2018): la conversión cm⁻¹ -> eV difiere en 8 partes por mil
# millones. Se compara en relativo con margen de sobra para eso y ninguno
# para un factor mal puesto.
REL_CODATA = 1e-6


def _ase_invcm():
    from ase.units import invcm
    return invcm


@pytest.mark.parametrize("nu", [
    [3800.0, 3700.0, 1600.0],                       # H2O
    [2100.0, 800.0, 600.0, 450.0, 300.0, 120.0],    # adsorbato
    [1000.0, 90.0, 45.0, 30.0, 500.0, 250.0],       # con modos blandos
])
def test_la_termoquimica_vibracional_concuerda_con_ase(nu):
    from ase.thermochemistry import HarmonicThermo

    from qekit.modules import thermochem as tc

    T = 298.15
    h = HarmonicThermo(np.array(nu) * _ase_invcm())
    zpe_ase = h.get_ZPE_correction()
    assert tc.zpe(nu) == pytest.approx(zpe_ase, rel=REL_CODATA)
    # se compara la energia interna ENTERA, no H_vib suelto: H_vib es la
    # diferencia de dos numeros del orden del ZPE y restarlos amplifica la
    # discrepancia de CODATA entre ASE y el paquete hasta hacerla visible.
    u = h.get_internal_energy(T, verbose=False)
    assert tc.zpe(nu) + tc.H_vib(nu, T) == pytest.approx(u, rel=REL_CODATA)
    assert tc.S_vib(nu, T) == pytest.approx(
        h.get_entropy(T, verbose=False), rel=REL_CODATA)


@pytest.mark.parametrize("nombre,simetria,geometria", [
    ("H2O", 2, "nonlinear"),
    ("CO2", 2, "linear"),
    ("NH3", 3, "nonlinear"),
    ("CH4", 12, "nonlinear"),
    ("N2", 2, "linear"),
    ("CO", 1, "linear"),
])
def test_la_entropia_de_gas_ideal_concuerda_con_ase(nombre, simetria, geometria):
    """Sackur-Tetrode + rotor rígido + vibraciones, molécula a molécula.

    El número de simetría entra como −k·ln(σ) y olvidarlo es el error
    clásico: el metano (σ = 12) se iría 0.06 eV a 300 K.
    """
    from ase.build import molecule
    from ase.thermochemistry import IdealGasThermo

    from qekit.modules import thermochem as tc

    at = molecule(nombre)
    nvib = 3 * len(at) - (5 if geometria == "linear" else 6)
    nu = np.linspace(600.0, 3200.0, nvib)
    T, p = 298.15, 101325.0
    ig = IdealGasThermo(vib_energies=nu * _ase_invcm(), geometry=geometria,
                        atoms=at, symmetrynumber=simetria, spin=0)
    nuestro = (tc.S_traslacional(float(at.get_masses().sum()), T, p)
               + tc.S_rotacional(at, T, simetria)
               + tc.S_vib(nu, T))
    assert nuestro == pytest.approx(ig.get_entropy(T, p, verbose=False),
                                    rel=REL_CODATA)


def test_los_momentos_de_inercia_son_los_de_ase():
    from ase.build import molecule

    from qekit.modules import thermochem as tc

    for nombre in ("H2O", "NH3", "CH4", "CO2"):
        at = molecule(nombre)
        assert tc.momentos_inercia(at) == pytest.approx(
            np.sort(at.get_moments_of_inertia()), abs=1e-9)


def test_el_numero_de_simetria_resta_k_ln_sigma():
    from ase.build import molecule

    from qekit.modules import thermochem as tc

    at = molecule("CH4")
    s1 = tc.S_rotacional(at, 300.0, 1)
    s12 = tc.S_rotacional(at, 300.0, 12)
    assert s1 - s12 == pytest.approx(tc.KB_EV * np.log(12.0), rel=1e-12)


# ----------------------------------------------------------------------
# Dinámica molecular contra límites analíticos
# ----------------------------------------------------------------------
def _tray(pos, celda, simbolos, dt=1.0):
    from qekit.modules import dynamics as dy
    return dy.Trayectoria(simbolos=list(simbolos), posiciones=pos,
                          celda=celda, dt=dt,
                          tiempos=np.arange(len(pos), dtype=float) * dt)


def test_la_gr_de_un_gas_ideal_vale_uno():
    """Si la normalización del cascarón está mal, esto no da 1."""
    from qekit.modules import dynamics as dy

    rng = np.random.default_rng(0)
    L, nat = 30.0, 900
    pos = rng.random((40, nat, 3)) * L
    tr = _tray(pos, np.eye(3) * L, ["Ar"] * nat)
    r, gr, _ = dy.rdf(tr, nbins=60)
    lejos = r > 3.0
    assert gr[lejos].mean() == pytest.approx(1.0, abs=0.02)
    assert gr[lejos].std() < 0.03


def test_el_numero_de_coordinacion_de_un_fcc_es_doce():
    from ase.build import bulk

    from qekit.modules import dynamics as dy

    rng = np.random.default_rng(1)
    at = bulk("Cu", "fcc", a=3.615, cubic=True).repeat(4)
    pos = at.get_positions()[None] + rng.normal(0, 0.05, (30, len(at), 3))
    tr = _tray(pos, at.cell.array, at.get_chemical_symbols())
    r, gr, _ = dy.rdf(tr, nbins=300)
    assert r[int(np.argmax(gr))] == pytest.approx(3.615 / np.sqrt(2), abs=0.03)
    _rmin, n = dy.coordinacion(r, gr, tr)
    assert n == pytest.approx(12.0, abs=0.4)


def test_la_difusion_recupera_el_D_de_un_movimiento_browniano():
    """MSD = 6·D·t. El 6 es el convenio tridimensional; con 2 o con 4 el
    número sale mal por un factor y nada lo delata."""
    from qekit.modules import dynamics as dy

    rng = np.random.default_rng(2)
    D_ref = 1.0e-5                       # cm²/s
    var = 6.0 * (D_ref * 10.0) / 3.0     # Å² por componente y por fs
    pos = np.cumsum(rng.normal(0.0, np.sqrt(var), (1200, 150, 3)), axis=0)
    tr = _tray(pos, np.eye(3) * 1e4, ["Li"] * 150)
    t, ms = dy.msd(tr)
    D, r2 = dy.difusion(t, ms)
    assert D == pytest.approx(D_ref, rel=0.06)
    assert r2 > 0.99
    # y la pendiente cruda, para fijar el 1/6 y el factor Å²/fs -> cm²/s
    pendiente = np.polyfit(t[100:500], ms[100:500], 1)[0]
    assert pendiente / 6.0 * 0.1 == pytest.approx(D_ref, rel=0.06)


def test_la_vdos_pone_el_pico_en_la_frecuencia_del_oscilador():
    """Fija el eje de frecuencias: 1/fs -> cm⁻¹ sin 2π de más ni de menos."""
    from qekit.modules import dynamics as dy

    rng = np.random.default_rng(3)
    nu = 500.0                                     # cm⁻¹
    f_fs = nu * 2.99792458e10 / 1e15
    t = np.arange(4000.0)
    osc = 0.1 * np.sin(2 * np.pi * f_fs * t)[:, None, None] * np.ones((1, 8, 3))
    pos = rng.random((1, 8, 3)) * 20.0 + osc
    tr = _tray(pos, np.eye(3) * 40.0, ["Si"] * 8)
    cm1, espectro = dy.vdos(tr)
    assert cm1[int(np.argmax(espectro))] == pytest.approx(nu, abs=5.0)


def test_desdoblar_deshace_exactamente_el_plegado():
    from qekit.modules import dynamics as dy

    rng = np.random.default_rng(4)
    celda = np.eye(3) * 10.0
    real = np.cumsum(rng.normal(0, 0.8, (300, 5, 3)), axis=0)
    plegado = real - np.round(real @ np.linalg.inv(celda)) @ celda
    des = dy.desdoblar(_tray(plegado, celda, ["H"] * 5))
    assert np.abs((des - des[0]) - (real - real[0])).max() < 1e-9


# ----------------------------------------------------------------------
# Transporte contra el desarrollo de Sommerfeld
# ----------------------------------------------------------------------
def _banda_modelo(p, mu, nk=20001, emax=20.0):
    """Estados con Σ(E) = w·v² ∝ E^p sobre una malla de energía uniforme."""
    from qekit.modules import transport as tp

    E = np.linspace(1e-6, emax, nk)
    run = tp.TransportRun()
    run.energies = E[:, None]
    run.velocities = np.zeros((nk, 1, 3))
    run.velocities[:, 0, 0] = np.sqrt(E ** p)
    run.weights = np.full(nk, 1.0 / nk)
    run.volume, run.fermi, run.nelec = 100.0, mu, 1.0
    return run


@pytest.mark.parametrize("p,T", [(0.5, 100.0), (1.5, 100.0), (3.0, 200.0)])
def test_el_seebeck_reproduce_la_formula_de_mott(p, T):
    """S = −(π²/3)(k_B²T/e)·dlnΣ/dE. Fija el signo y el prefactor a la vez:
    con el signo al revés un tipo p daría Seebeck negativo."""
    from qekit.modules import transport as tp

    mu = 8.0
    run = tp.compute(_banda_modelo(p, mu), T=[T], mu=[mu])
    mott = -(np.pi ** 2 / 3.0) * tp.KB_EV ** 2 * T * p / mu
    assert run.seebeck[0, 0, 0, 0] == pytest.approx(mott, rel=3e-3)


@pytest.mark.parametrize("p,T", [(0.5, 100.0), (3.0, 200.0)])
def test_el_numero_de_lorenz_tiende_a_L0_en_el_limite_degenerado(p, T):
    """κ_e/(σT) -> L₀ sea cual sea la forma de la banda. Si a κ_e le faltara
    la resta de S²σT, esto no daría L₀."""
    from qekit.modules import transport as tp

    mu = 8.0
    run = tp.compute(_banda_modelo(p, mu), T=[T], mu=[mu])
    L = run.kappa_e[0, 0, 0, 0] / (run.sigma[0, 0, 0, 0] * T)
    assert L == pytest.approx(tp.L0_SOMMERFELD, rel=3e-3)


def test_L0_es_el_exacto_y_no_el_2_44e_8_de_los_libros():
    from qekit.modules import transport as tp

    assert tp.L0_SOMMERFELD == pytest.approx(
        (np.pi ** 2 / 3.0) * tp.KB_EV ** 2, rel=1e-15)
    assert abs(tp.L0_SOMMERFELD - 2.44e-8) / tp.L0_SOMMERFELD > 1e-3


# ----------------------------------------------------------------------
# Ópticas: identidades exactas
# ----------------------------------------------------------------------
def _lorentz(E, E0=4.0, f=12.0, g=0.6):
    den = (E0 ** 2 - E ** 2) ** 2 + (g * E) ** 2
    return 1.0 + f * (E0 ** 2 - E ** 2) / den, f * g * E / den


def test_n_y_k_son_la_raiz_compleja_de_epsilon():
    from qekit.modules import optics as op

    E = np.linspace(0.02, 40.0, 4000)
    e1, e2 = _lorentz(E)
    run = op.OpticsRun()
    run.energies, run.eps1, run.eps2 = E, e1, e2
    d = op.derived(run)
    n, k = d["n"], d["k"]
    assert np.abs(n ** 2 - k ** 2 - e1).max() < 1e-11
    assert np.abs(2 * n * k - e2).max() < 1e-11
    # Fresnel a incidencia normal
    nc = n + 1j * k
    assert np.abs(d["R"] - np.abs((nc - 1) / (nc + 1)) ** 2).max() < 1e-12
    # α = 4πk/λ, con λ sacada de ħc
    lam_cm = 2 * np.pi * op.HBAR_C_EV_CM / E
    assert np.abs(d["alpha"] - 4 * np.pi * k / lam_cm).max() / d["alpha"].max() < 1e-12


def test_kramers_kronig_reconstruye_un_oscilador_de_lorentz():
    from qekit.modules import optics as op

    E = np.linspace(0.02, 40.0, 8000)
    e1, e2 = _lorentz(E)
    e1_kk = op.kramers_kronig(E, e2)
    dentro = (E > 0.5) & (E < 30.0)
    assert np.abs(e1_kk - e1)[dentro].max() < 0.03
    assert e1_kk[-1] == pytest.approx(1.0, abs=0.02)


# ----------------------------------------------------------------------
# Termoelásticas derivadas
# ----------------------------------------------------------------------
def test_la_kappa_de_slack_no_depende_de_que_celda_se_le_de():
    """El mismo cristal, descrito con la primitiva o con la convencional,
    tiene que dar la misma κ. La n del modelo cuenta ramas ópticas, así que
    es la de la celda PRIMITIVA; usar la del fichero daba 95.3 W/mK con la
    primitiva del silicio y 37.8 con la convencional."""
    from qekit.modules import derived

    B, G = 97.9, 66.5                       # GPa, silicio
    prim = derived.analyze(B, G, [28.0855] * 2, 40.04, natoms=2,
                           n_primitiva=2, T=300.0)
    conv = derived.analyze(B, G, [28.0855] * 8, 160.16, natoms=8,
                           n_primitiva=2, T=300.0)
    assert conv.kappa_slack == pytest.approx(prim.kappa_slack, rel=1e-12)
    # lo que sí es cell-independent por construcción
    assert conv.rho == pytest.approx(prim.rho, rel=1e-12)
    assert conv.theta_D == pytest.approx(prim.theta_D, rel=1e-12)


def test_la_debye_elastica_del_silicio_sale_donde_la_mide_el_experimento():
    """Con B y G experimentales, θ_D = 645 K. Ata de una vez las tres
    fórmulas encadenadas: densidad, velocidades de Anderson y (6π²n)^(1/3)."""
    from qekit.modules import derived

    r = derived.analyze(97.9, 66.5, [28.0855] * 2, 40.04, natoms=2)
    assert r.rho == pytest.approx(2329.0, abs=5.0)
    assert r.theta_D == pytest.approx(645.0, abs=15.0)


def test_las_relaciones_elasticas_derivadas_son_identidades():
    from qekit.modules import derived

    B, G, rho = 97.9, 66.5, 2329.0
    v_l, v_t, v_m = derived.sound_velocities(B, G, rho)
    assert v_l == pytest.approx(np.sqrt((B + 4 * G / 3) * 1e9 / rho), rel=1e-12)
    assert v_t == pytest.approx(np.sqrt(G * 1e9 / rho), rel=1e-12)
    assert 3.0 / v_m ** 3 == pytest.approx(2 / v_t ** 3 + 1 / v_l ** 3, rel=1e-12)
    nu = derived.poisson_ratio(B, G)
    # ν y E salen del mismo par (B, G): se comprueba contra E = 9BG/(3B+G)
    E = 9 * B * G / (3 * B + G)
    assert nu == pytest.approx(E / (2 * G) - 1.0, rel=1e-12)
    # Belomestnykh-Tesleva
    assert derived.gruneisen_from_poisson(nu) == pytest.approx(
        1.5 * (1 + nu) / (2 - 3 * nu), rel=1e-12)


def test_la_debye_por_momento_de_la_dos_usa_el_segundo_momento():
    """θ_D = (ħ/k_B)·√(5/3·⟨ω²⟩) sobre una DOS de Debye analítica, donde
    ⟨ω²⟩ = 3/5·ω_max² y por tanto θ_D sale exactamente ħω_max/k_B."""
    from qekit.modules import derived

    w = np.linspace(1.0, 600.0, 4000)          # cm⁻¹
    dos = w ** 2                               # DOS de Debye
    theta = derived.debye_from_dos(w, dos, natoms=1)
    w_rad = 600.0 * 2 * np.pi * 2.99792458e10
    esperado = derived.HBAR / derived.KB * w_rad
    assert theta == pytest.approx(esperado, rel=2e-3)
    # y la normalización no puede cambiar el resultado
    assert derived.debye_from_dos(w, 7.5 * dos, natoms=1) == pytest.approx(
        theta, rel=1e-12)
    assert trapezoid(dos, w) > 0


# ----------------------------------------------------------------------
# Superconductividad: Allen-Dynes contra su caso de calibración
# ----------------------------------------------------------------------
def test_allen_dynes_reproduce_los_7_2_K_del_plomo():
    """El plomo (λ = 1.55, ω_log = 56 K, μ* = 0.10) es el caso con el que se
    calibró la fórmula. Ata a la vez el 1.04, el 0.62, el 1.2 y el factor f1:
    sin f1 la fórmula desnuda da 6.6 K."""
    from qekit.modules.elph import allen_dynes

    assert allen_dynes(1.55, 56.0, 0.10) == pytest.approx(7.2, abs=0.1)
    assert allen_dynes(1.55, 56.0, 0.10, correcciones=False) == pytest.approx(
        6.6, abs=0.1)


def test_el_informe_avisa_cuando_lambda_sale_del_ajuste():
    from qekit.modules import elph

    run = elph.ElPhRun()
    run.sigmas = np.array([0.02])
    run.omega_log = np.array([400.0])
    run.mustar = 0.10
    run.lambdas = np.array([2.4])
    assert "beyond where Allen-Dynes was fitted" in elph.report(run)
    run.lambdas = np.array([0.8])
    assert "beyond where Allen-Dynes was fitted" not in elph.report(run)


# ----------------------------------------------------------------------
# Difracción: hasta dónde vale la inversión de Mott-Bethe
# ----------------------------------------------------------------------
def test_los_factores_de_dispersion_son_mott_bethe_y_no_cromer_mann():
    """El 41.78214 es 8π²a₀. Concuerda con Cromer-Mann por debajo de
    s ≈ 0.5 Å⁻¹ y se separa después: es una condición de validez, no un
    error, pero hay que saberla."""
    from qekit.modules import xrd

    assert 8 * np.pi ** 2 * 0.529177210903 == pytest.approx(41.78214, abs=1e-4)
    # Cromer-Mann del silicio (International Tables C, tabla 6.1.1.4)
    a = np.array([6.2915, 3.0353, 1.9891, 1.5410])
    b = np.array([2.4386, 32.3337, 0.6785, 81.6937])
    c = 1.1407
    ab = xrd.scattering_params("Si")
    for s, tol in [(0.0, 0.01), (0.2, 0.02), (0.5, 0.05), (1.5, None)]:
        s2 = s ** 2
        mb = 14 - 41.78214 * s2 * np.sum(ab[:, 0] * np.exp(-ab[:, 1] * s2))
        cm = float(np.sum(a * np.exp(-b * s2)) + c)
        if tol is None:                      # a s alta se separan de verdad
            assert abs(mb - cm) / cm > 0.15
        else:
            assert abs(mb - cm) / cm < tol


def test_el_difractograma_avisa_solo_cuando_llega_a_s_alta():
    from ase.build import bulk

    from qekit.modules import xrd

    at = bulk("Si", "diamond", a=5.431, cubic=True)
    rango = (5.0, 140.0)
    cu = xrd.report(xrd.compute(at, wavelength="CuKa", two_theta_range=rango))
    ag = xrd.report(xrd.compute(at, wavelength="AgKa", two_theta_range=rango))
    assert "Mott-Bethe" not in cu          # con Cu no se llega: s <= 1/lambda
    assert "Mott-Bethe" in ag
    # y la nota de Scherrer solo sale si se pidió un tamaño de cristalito
    pat = xrd.broaden(xrd.compute(at, two_theta_range=rango),
                      two_theta_range=rango, size_nm=20.0)
    assert "Scherrer with K = 0.9" in xrd.report(pat)


# ----------------------------------------------------------------------
# Termodinámica armónica: phonopy, consistencia interna y entre módulos
# ----------------------------------------------------------------------
#: phonopy deriva su k_B de su propio juego de unidades y sale 5.8e-7 por
#: encima del exacto de CODATA 2018 (k_B/e con las dos constantes que
#: definen el SI). La comparación se hace con ese margen, no más.
REL_PHONOPY = 5e-6

MODOS_CM1 = np.array([1000.0, 620.0, 450.0, 310.0, 180.0, 95.0])


def test_las_constantes_de_conversion_son_las_mismas_en_todo_el_paquete():
    """Tres módulos convierten cm⁻¹ a eV por su cuenta. Si uno se desvía,
    la energía libre de `qha` y la de `thermochem` dejan de ser la misma
    magnitud sin que nada falle."""
    from qekit.modules import phonons as ph
    from qekit.modules import qha
    from qekit.modules import thermochem as tc

    exacto = 299792458 * 100 * 4.135667696e-15      # c[cm/s] * h[eV s]
    for valor in (qha.CM1_EV, ph.CM1_TO_EV, float(tc.cm1_a_eV(1.0))):
        assert valor == pytest.approx(exacto, rel=1e-10)
    # k_B en eV/K es exacto en el SI: k_B/e con las dos constantes que lo definen
    assert tc.KB_EV == pytest.approx(1.380649e-23 / 1.602176634e-19, rel=1e-10)


@pytest.mark.parametrize("T", [50.0, 100.0, 300.0, 800.0])
def test_la_capacidad_calorifica_armonica_concuerda_con_phonopy(T):
    from phonopy.phonon.thermal_properties import mode_cv

    from qekit.modules import qha
    from qekit.modules import thermochem as tc

    ref = float(np.sum(mode_cv(T, MODOS_CM1 * qha.CM1_EV)))     # eV/K
    assert tc.Cv_vib(MODOS_CM1, T) == pytest.approx(ref, rel=REL_PHONOPY)
    # cv_modos va en meV/K
    assert qha.cv_modos(MODOS_CM1, T) / 1000.0 == pytest.approx(
        ref, rel=REL_PHONOPY)


@pytest.mark.parametrize("T", [50.0, 300.0, 800.0])
def test_qha_y_thermochem_dan_la_misma_energia_libre(T):
    """Dos implementaciones independientes dentro del paquete. Si divergen,
    una de las dos está mal y nada lo delataría."""
    from qekit.modules import qha
    from qekit.modules import thermochem as tc

    esperado = (tc.zpe(MODOS_CM1) + tc.H_vib(MODOS_CM1, T)
                - T * tc.S_vib(MODOS_CM1, T))
    assert qha.f_vib(MODOS_CM1, T) == pytest.approx(esperado, rel=1e-9)


@pytest.mark.parametrize("T", [100.0, 300.0, 800.0])
def test_las_funciones_termodinamicas_son_derivadas_unas_de_otras(T):
    """C_v = dU/dT y S = −dF/dT, por diferencias finitas. Es lo que atrapa
    un factor perdido en una sola de las tres."""
    from qekit.modules import qha
    from qekit.modules import thermochem as tc

    h = 1e-3 * T
    dU = ((tc.H_vib(MODOS_CM1, T + h) - tc.H_vib(MODOS_CM1, T - h)) / (2 * h))
    assert dU == pytest.approx(tc.Cv_vib(MODOS_CM1, T), rel=1e-5)
    dF = ((qha.f_vib(MODOS_CM1, T + h) - qha.f_vib(MODOS_CM1, T - h)) / (2 * h))
    assert -dF == pytest.approx(tc.S_vib(MODOS_CM1, T), rel=1e-5)


def test_dulong_petit():
    from qekit.modules import thermochem as tc

    assert tc.Cv_vib(MODOS_CM1, 40000.0) == pytest.approx(
        len(MODOS_CM1) * tc.KB_EV, rel=1e-3)


def test_la_termodinamica_por_DOS_reproduce_la_suma_sobre_modos():
    """`phonons.thermodynamics` integra una DOS; `qha` suma modos. Sobre una
    DOS que son esos mismos modos ensanchados, tienen que coincidir."""
    from qekit.modules import phonons as ph
    from qekit.modules import qha
    from qekit.modules import thermochem as tc

    w = np.linspace(1.0, 1400.0, 40000)
    sigma = 2.0
    g = np.zeros_like(w)
    for w0 in MODOS_CM1:
        g += np.exp(-0.5 * ((w - w0) / sigma) ** 2) / (sigma * np.sqrt(2 * np.pi))
    run = ph.PhononRun()
    run.dos_w, run.dos = w, g
    Ts = np.array([100.0, 300.0, 800.0])
    res = ph.thermodynamics(run, T=Ts, natoms=2)          # 3N = 6 modos
    assert res["ZPE"] == pytest.approx(tc.zpe(MODOS_CM1), rel=1e-3)
    for i, T in enumerate(Ts):
        assert res["F"][i] == pytest.approx(qha.f_vib(MODOS_CM1, T), rel=3e-3)
        assert res["Cv"][i] == pytest.approx(tc.Cv_vib(MODOS_CM1, T), rel=3e-3)
        # la relación que define S en ese módulo
        assert res["F"][i] == pytest.approx(res["U"][i] - T * res["S"][i],
                                            rel=1e-12)


# ----------------------------------------------------------------------
# Masa efectiva, ecuación de estado, elasticidad y Tauc
# ----------------------------------------------------------------------
def test_la_masa_efectiva_de_un_electron_libre_es_uno():
    """E = ħ²k²/2mₑ da a = 3.80998212 eV·Å². Si el 2 de m* = ħ²/(2a) se
    perdiera, toda masa efectiva del paquete saldría con un factor 2."""
    from qekit.modules import effmass

    assert effmass._mass_from_quadratic(3.80998212) == pytest.approx(1.0, rel=1e-7)
    # y la constante es la identidad, no un literal de una edición de CODATA
    assert effmass.HBAR2_OVER_ME == pytest.approx(
        27.211386245988 * 0.529177210903 ** 2, rel=1e-8)
    # el signo distingue hueco de electrón
    assert effmass._mass_from_quadratic(-3.80998212) < 0


@pytest.mark.parametrize("ecuacion,funcion", [
    ("birch-murnaghan", "birch_murnaghan"),
    ("murnaghan", "murnaghan"),
    ("vinet", "vinet"),
])
def test_las_ecuaciones_de_estado_recuperan_sus_propios_parametros(ecuacion, funcion):
    """Se genera la curva con (V0, B0, B0') conocidos y se vuelve a ajustar.
    Cierra el lazo entre la fórmula, el ajuste y la conversión eV/Å³ -> GPa."""
    from qekit.modules import eos

    V0, E0, B0_GPa, Bp = 40.0, -10.0, 90.0, 4.2
    EV_A3_GPA = 160.21766208
    V = np.linspace(0.90 * V0, 1.10 * V0, 11)
    E = getattr(eos, funcion)(V, E0, V0, B0_GPa / EV_A3_GPA, Bp)
    run = eos.EOSRun()
    run.volumes, run.energies = V, E
    f = eos.fit(run, equation=ecuacion)
    assert f.V0 == pytest.approx(V0, rel=1e-6)
    assert f.B0 == pytest.approx(B0_GPa, rel=1e-5)
    assert f.Bp == pytest.approx(Bp, rel=1e-5)


def test_birch_murnaghan_concuerda_con_ase():
    from ase.eos import EquationOfState

    from qekit.modules import eos

    V0, E0, B0_GPa, Bp = 40.0, -10.0, 90.0, 4.2
    EV_A3_GPA = 160.21766208
    V = np.linspace(0.90 * V0, 1.10 * V0, 11)
    E = eos.birch_murnaghan(V, E0, V0, B0_GPa / EV_A3_GPA, Bp)
    v0, _e0, b0 = EquationOfState(list(V), list(E), eos="birchmurnaghan").fit()
    assert v0 == pytest.approx(V0, rel=1e-6)
    assert b0 * EV_A3_GPA == pytest.approx(B0_GPa, rel=1e-4)


def test_en_un_solido_isotropo_voigt_y_reuss_coinciden():
    """Con C construido de λ y μ, las cotas tienen que cerrarse: si no,
    hay un error en la inversión de C o en los índices de Voigt."""
    from qekit.modules import elastic

    lam, mu = 60.0, 40.0                              # GPa
    C = np.zeros((6, 6))
    C[:3, :3] = lam
    C[0, 0] = C[1, 1] = C[2, 2] = lam + 2 * mu
    C[3, 3] = C[4, 4] = C[5, 5] = mu
    m = elastic.moduli(C)
    B, G = lam + 2 * mu / 3.0, mu
    assert m.B_voigt == pytest.approx(B, rel=1e-12)
    assert m.B_reuss == pytest.approx(B, rel=1e-12)
    assert m.G_voigt == pytest.approx(G, rel=1e-12)
    assert m.G_reuss == pytest.approx(G, rel=1e-12)
    assert m.E == pytest.approx(9 * B * G / (3 * B + G), rel=1e-12)
    assert m.nu == pytest.approx((3 * B - 2 * G) / (2 * (3 * B + G)), rel=1e-12)
    assert abs(m.anisotropy) < 1e-9        # A^U = 0 exactamente si es isótropo


@pytest.mark.parametrize("clase,potencia", [("direct", 0.5), ("indirect", 2.0)])
def test_tauc_recupera_un_gap_puesto_a_mano(clase, potencia):
    """ε₂ construida para que (αhν)^(1/r) sea exactamente una recta que corta
    en Eg. Fija el exponente de cada clase y el corte con el eje."""
    from qekit.modules import optics as op

    Eg = 1.85
    E = np.arange(0.05, 12.0, 0.005)
    borde = np.maximum(E - Eg, 0.0)
    e2 = np.where(E > Eg, borde ** potencia / E ** 2, 0.0)
    e2 = e2 * (30.0 if clase == "direct" else 0.6)
    run = op.OpticsRun()
    run.energies, run.eps1, run.eps2 = E, np.full_like(E, 12.0), e2
    gap, pendiente, _ventana, _y = op.tauc_gap(run, kind=clase)
    assert gap == pytest.approx(Eg, abs=0.03)
    assert pendiente > 0


# ----------------------------------------------------------------------
# Defectos cargados, conductancia y casco convexo
# ----------------------------------------------------------------------
@pytest.mark.parametrize("red,alpha_ref", [
    ("sc", 2.8373), ("bcc", 2.8883), ("fcc", 2.8883),
])
def test_las_constantes_de_madelung_son_las_de_makov_payne(red, alpha_ref):
    """α con L = V^(1/3) para una carga puntual en fondo neutralizante.
    Son los tres números de la tabla de Makov-Payne y no dependen de nada
    que haya que ir a buscar: la suma de Ewald los produce."""
    from qekit.modules import defects

    a = 5.0
    celdas = {
        "sc": np.eye(3) * a,
        "bcc": a * np.array([[-.5, .5, .5], [.5, -.5, .5], [.5, .5, -.5]]),
        "fcc": a * np.array([[0, .5, .5], [.5, 0, .5], [.5, .5, 0]]),
    }
    assert defects.constante_madelung(celdas[red]) == pytest.approx(
        alpha_ref, abs=5e-4)


def test_la_constante_de_madelung_no_depende_del_tamano_de_la_celda():
    from qekit.modules import defects

    a5 = defects.constante_madelung(np.eye(3) * 5.0)
    a17 = defects.constante_madelung(np.eye(3) * 17.0)
    assert a5 == pytest.approx(a17, rel=1e-7)


def test_la_correccion_de_imagen_escala_como_q2_sobre_epsilon_L():
    from qekit.modules import defects

    c = np.eye(3) * 10.0
    def mp(q, celda, eps):
        return defects.correccion_imagen(q, celda, eps, "makov-payne")["E_mp"]
    base = mp(1, c, 10.0)
    assert mp(2, c, 10.0) / base == pytest.approx(4.0, rel=1e-9)
    assert mp(1, np.eye(3) * 20.0, 10.0) / base == pytest.approx(0.5, rel=1e-6)
    assert mp(1, c, 20.0) / base == pytest.approx(0.5, rel=1e-12)
    # e²/(4πε₀) en eV·Å
    assert defects.KE == pytest.approx(14.399645, rel=1e-7)


def test_el_cuanto_de_conductancia_es_2e2_sobre_h():
    from qekit.modules import ballistic

    exacto = 2 * (1.602176634e-19) ** 2 / 6.62607015e-34
    assert ballistic.G0 == pytest.approx(exacto, rel=1e-9)
    assert ballistic.R0 == pytest.approx(1.0 / exacto, rel=1e-9)
    assert ballistic.R0 == pytest.approx(12906.4, abs=0.1)


def test_el_casco_convexo_binario_da_las_distancias_exactas():
    """Cinco fases con la respuesta calculable a mano: A3B queda 0.1 eV/át
    por encima del segmento A–AB, y AB3 0.2 por encima del AB–B."""
    from qekit.modules import thermo

    filas = [("A", {"A": 1}, 0.0), ("B", {"B": 1}, 0.0),
             ("A3B", {"A": 3, "B": 1}, -0.4 * 4),
             ("AB", {"A": 1, "B": 1}, -1.0 * 2),
             ("AB3", {"A": 1, "B": 3}, -0.3 * 4)]
    res = thermo.from_table(filas, elementos=["A", "B"])
    esperado = {"A": 0.0, "B": 0.0, "A3B": 0.1, "AB": 0.0, "AB3": 0.2}
    for f in res.fases:
        assert f.e_hull == pytest.approx(esperado[f.nombre], abs=1e-9)
        assert f.en_casco == (esperado[f.nombre] < 1e-9)


def test_el_casco_ternario_interpola_en_el_simplex_correcto():
    """A2BC está sobre la línea A–ABC: 0.75 de ABC y 0.25 de A dan −1.125,
    así que su distancia al casco es exactamente 0.225 eV/átomo."""
    from qekit.modules import thermo

    filas = [("A", {"A": 1}, 0.0), ("B", {"B": 1}, 0.0), ("C", {"C": 1}, 0.0),
             ("ABC", {"A": 1, "B": 1, "C": 1}, -1.5 * 3),
             ("A2BC", {"A": 2, "B": 1, "C": 1}, -0.9 * 4)]
    res = thermo.from_table(filas, elementos=["A", "B", "C"])
    por_nombre = {f.nombre: f for f in res.fases}
    assert por_nombre["ABC"].e_hull == pytest.approx(0.0, abs=1e-9)
    assert por_nombre["ABC"].en_casco
    assert por_nombre["A2BC"].e_hull == pytest.approx(0.225, abs=1e-9)
    assert not por_nombre["A2BC"].en_casco
