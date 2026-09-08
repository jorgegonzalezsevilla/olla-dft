"""Número de Lorenz y transporte por canal de espín."""

import numpy as np
import pytest

from qekit.modules import transport as tr

HBAR = 6.582119569e-16                       # eV s
ME = 0.51099895e6 / (2.99792458e8) ** 2      # eV s^2 / m^2


def gas_de_electrones(n=20, a=4.0, T=300.0):
    """Gas de electrones libres: E = ħ²k²/2m, v = ħk/m.

    Para un metal degenerado el numero de Lorenz vale exactamente
    L0 = (pi^2/3)(k_B/e)^2 = 2.44e-8. Es el unico caso con respuesta
    analitica, asi que es EL test del modulo.
    """
    cell = np.eye(3) * a
    recip = 2 * np.pi * np.linalg.inv(cell).T
    ks = (np.arange(n) + 0.5) / n - 0.5
    K = np.array(np.meshgrid(ks, ks, ks, indexing="ij")).reshape(3, -1).T
    kSI = (K @ recip) * 1e10
    E = HBAR ** 2 * np.sum(kSI ** 2, axis=1) / (2 * ME)
    v = HBAR * kSI / ME
    run = tr.TransportRun(volume=a ** 3, nelec=2.0,
                          fermi=float(np.median(E)), grid=(n, n, n))
    run.energies = E[:, None]
    run.velocities = v[:, None, :]
    run.weights = np.full(len(E), 2.0 / len(E))
    return tr.compute(run, T=[T], mu=np.array([run.fermi]))


def test_el_gas_de_electrones_libres_da_L0():
    """La comprobación decisiva: L -> 2.44e-8 en un metal degenerado."""
    run = gas_de_electrones()
    L = float(tr.lorenz(run, 0)[0])
    assert L / tr.L0_SOMMERFELD == pytest.approx(1.0, abs=0.10), \
        f"L/L0 = {L / tr.L0_SOMMERFELD:.3f}"


def test_el_seebeck_de_un_metal_es_pequeno():
    run = gas_de_electrones()
    S = float(np.trace(run.seebeck[0][0]) / 3.0)
    assert abs(S) < 60e-6, "un metal degenerado da decenas de uV/K, no cientos"


def test_L_no_depende_de_tau():
    """Multiplicar sigma y kappa por el mismo factor no cambia L."""
    run = gas_de_electrones()
    L1 = float(tr.lorenz(run, 0)[0])
    run.sigma = run.sigma * 7.3
    run.kappa_e = run.kappa_e * 7.3
    assert float(tr.lorenz(run, 0)[0]) == pytest.approx(L1, rel=1e-12)


def test_L_no_depende_de_la_temperatura_en_un_metal():
    a = float(tr.lorenz(gas_de_electrones(T=200.0), 0)[0])
    b = float(tr.lorenz(gas_de_electrones(T=500.0), 0)[0])
    assert a == pytest.approx(b, rel=0.10)


def test_el_metal_no_dispara_la_alarma_de_cancelacion():
    run = gas_de_electrones()
    assert float(tr.cancelacion(run, 0)[0]) > 0.1
    assert "DO NOT TRUST" not in tr.report_lorenz(run)


# ----------------------------------------------------------------------
# cancelación catastrófica
# ----------------------------------------------------------------------
def _run_falso(sig, See, kap, T=300.0):
    r = tr.TransportRun(volume=40.0, fermi=0.0, grid=(4, 4, 4))
    n = len(sig)
    r.mu = np.linspace(-1, 1, n)
    r.T = np.array([T])
    eye = np.eye(3)
    r.sigma = np.array([[s * eye for s in sig]])
    r.seebeck = np.array([[s * eye for s in See]])
    r.kappa_e = np.array([[k * eye for k in kap]])
    r.carriers = np.zeros((1, n))
    return r


def test_detecta_la_cancelacion_catastrofica():
    """kappa_e = kappa0 - S^2 sigma T; si de la resta queda el 0.003 %, es ruido."""
    sig, S, T = 6.14e15, 976.5e-6, 300.0
    k0 = S ** 2 * sig * T
    kap = k0 * 3.4e-5                       # lo que sobrevive de verdad
    r = _run_falso([sig] * 3, [S] * 3, [kap] * 3)
    c = tr.cancelacion(r, 0)
    assert c[1] < 1e-3
    assert "DO NOT TRUST" in tr.report_lorenz(r)


def test_sin_cancelacion_no_avisa():
    sig, S, T = 1e20, 20e-6, 300.0
    kap = tr.L0_SOMMERFELD * sig * T
    r = _run_falso([sig] * 3, [S] * 3, [kap] * 3)
    txt = tr.report_lorenz(r)
    assert "DO NOT TRUST" not in txt
    assert "Wiedemann-Franz holds within 15" in txt


def test_por_encima_de_L0_habla_de_bipolar():
    sig, S, T = 1e20, 20e-6, 300.0
    r = _run_falso([sig] * 3, [S] * 3, [3.0 * tr.L0_SOMMERFELD * sig * T] * 3)
    assert "BIPOLAR" in tr.report_lorenz(r)


def test_por_debajo_de_L0_habla_de_no_degenerado():
    sig, S, T = 1e20, 20e-6, 300.0
    r = _run_falso([sig] * 3, [S] * 3, [0.6 * tr.L0_SOMMERFELD * sig * T] * 3)
    txt = tr.report_lorenz(r)
    assert "is NOT degenerate" in txt and "0.76" in txt


# ----------------------------------------------------------------------
# dos canales de espín
# ----------------------------------------------------------------------
def _espin(sig_up, sig_dw, S_up, S_dw):
    up = _run_falso([sig_up] * 3, [S_up] * 3, [1.0] * 3)
    dw = _run_falso([sig_dw] * 3, [S_dw] * 3, [1.0] * 3)
    up.fermi = dw.fermi = 0.0
    return tr.TransporteEspin(up=up, dw=dw, it=0)


def test_las_conductancias_se_suman():
    te = _espin(3.0, 1.0, 10e-6, 10e-6)
    assert te.sigma_total[1] == pytest.approx(4.0)


def test_el_seebeck_se_pesa_por_la_conductancia():
    """Media pesada, no aritmetica: 3:1 de conductancia da (3*100+1*(-100))/4."""
    te = _espin(3.0, 1.0, 100e-6, -100e-6)
    assert te.seebeck_total[1] * 1e6 == pytest.approx(50.0)
    media_mala = 0.0
    assert te.seebeck_total[1] * 1e6 != pytest.approx(media_mala)


def test_un_canal_que_no_conduce_no_aporta_termopotencia():
    te = _espin(1.0, 0.0, 10e-6, 9999e-6)
    assert te.seebeck_total[1] * 1e6 == pytest.approx(10.0)


def test_la_polarizacion_va_de_menos_uno_a_uno():
    assert _espin(1.0, 0.0, 0, 0).polarizacion[1] == pytest.approx(1.0)
    assert _espin(0.0, 1.0, 0, 0).polarizacion[1] == pytest.approx(-1.0)
    assert _espin(1.0, 1.0, 0, 0).polarizacion[1] == pytest.approx(0.0)


def test_un_medio_metal_se_reconoce():
    te = _espin(1.0, 0.001, 20e-6, 20e-6)
    assert "half-metal" in tr.report_espin(te)


def test_sin_polarizacion_lo_dice():
    te = _espin(1.0, 1.0, 20e-6, 20e-6)
    assert "adds nothing" in tr.report_espin(te)


def test_la_termopotencia_de_espin_es_la_diferencia():
    te = _espin(1.0, 1.0, 80e-6, 20e-6)
    assert te.seebeck_de_espin[1] * 1e6 == pytest.approx(60.0)


# ----------------------------------------------------------------------
# Regresiones: la transformación de las velocidades y el convenio de pesos
# ----------------------------------------------------------------------
def test_las_velocidades_son_correctas_en_una_celda_hexagonal():
    """La red recíproca hexagonal NO es simétrica, y ahí se veía el fallo.

    Con E = kx² la velocidad tiene que ir toda en x. Aplicar B⁻¹ en vez de
    B⁻ᵀ metía una v_y = −v_x/2 espuria y subía |v|² un 25 %. En cúbica no se
    notaba porque ahí B sí es simétrica, que es por lo que el resto de las
    pruebas del módulo no lo cazaban.
    """
    a, c, n = 2.46, 10.0, 12
    cell = np.array([[a, 0.0, 0.0],
                     [-a / 2, a * np.sqrt(3) / 2, 0.0],
                     [0.0, 0.0, c]])
    recip = 2 * np.pi * np.linalg.inv(cell).T
    f = (np.arange(n) + 0.5) / n - 0.5
    F = np.array(np.meshgrid(f, f, f, indexing="ij")).reshape(3, -1).T
    kcart = F @ recip                                   # Å⁻¹
    E = (kcart[:, 0] ** 2).reshape(n, n, n, 1)          # eV, E = kx²

    v = tr._fd_derivative(E, cell, (n, n, n))            # (n,n,n,1,3)
    esperado = 2.0 * kcart[:, 0].reshape(n, n, n) * (tr.ANG_M / tr.HBAR_EVS)

    assert np.allclose(v[..., 1], 0.0), \
        "aparece una componente v_y que no existe: falta la traspuesta"
    assert np.allclose(v[..., 2], 0.0)
    # El módulo se compara en el INTERIOR del eje f1: ahí la diferencia
    # central es exacta para una parábola. En la primera y la última capa el
    # np.pad(mode="wrap") envuelve una función que no es periódica (una banda
    # de verdad sí lo es), y eso es artefacto de la prueba, no del módulo.
    assert np.allclose(v[1:-1, :, :, 0, 0], esperado[1:-1], rtol=1e-9)


def test_las_velocidades_no_cambian_al_girar_la_celda():
    """Invariancia física: girar cristal y malla no cambia el módulo de v."""
    n = 8
    cell = np.array([[3.0, 0.0, 0.0], [1.0, 3.5, 0.0], [0.0, 0.4, 4.0]])
    rng = np.random.default_rng(0)
    E = rng.normal(size=(n, n, n, 2))
    th = 0.7
    R = np.array([[np.cos(th), -np.sin(th), 0.0],
                  [np.sin(th), np.cos(th), 0.0],
                  [0.0, 0.0, 1.0]])

    v1 = tr._fd_derivative(E, cell, (n, n, n))
    v2 = tr._fd_derivative(E, cell @ R.T, (n, n, n))
    assert np.allclose(np.linalg.norm(v1, axis=-1),
                       np.linalg.norm(v2, axis=-1))


class _XMLFalso:
    """Lo mínimo que `load` mira de un nscf, para no depender de datos de QE."""

    def __init__(self, nspin, n=4, a=4.0):
        f = np.arange(n) / n
        self.cell = np.eye(3) * a
        self.calculation = "nscf"
        self.volume = a ** 3
        self.nelec = 2.0
        self.fermi = 0.0
        self.kpoints_frac = np.array(
            np.meshgrid(f, f, f, indexing="ij")).reshape(3, -1).T
        nk = len(self.kpoints_frac)
        self.eigenvalues = np.zeros((nspin, nk, 1))


def test_load_pone_la_degeneracion_de_espin_en_los_pesos(monkeypatch):
    """Los pesos siguen el convenio de QE: suman 2 sin polarizar, 1 por canal.

    Normalizarlos siempre a 1 dejaba σ/τ y κ_e/τ a la MITAD en el caso sin
    polarizar, que es el de por omisión.
    """
    for nspin, total in ((1, 2.0), (2, 1.0)):
        monkeypatch.setattr(tr.qeout, "read_xml",
                            lambda _p, _n=nspin: _XMLFalso(_n))
        run = tr.load("da igual.xml")
        assert run.nspin == nspin
        assert float(run.weights.sum()) == pytest.approx(total), \
            f"con nspin={nspin} los pesos tienen que sumar {total}"


def test_sigma_lleva_la_degeneracion_de_espin():
    """σ/τ es LINEAL en los pesos, y los pesos llevan el espín.

    Sin polarizar cada estado aloja dos electrones (los pesos suman 2). Si
    se normalizaran a 1, σ/τ y κ_e/τ saldrían a la mitad; S y el número de
    Lorenz no lo notarían por ser cocientes, que es justo lo que escondía
    el fallo.
    """
    run = gas_de_electrones()
    mitad = tr.TransportRun(volume=run.volume, nelec=run.nelec,
                            fermi=run.fermi, grid=run.grid)
    mitad.energies, mitad.velocities = run.energies, run.velocities
    mitad.weights = run.weights / 2.0
    mitad = tr.compute(mitad, T=[300.0], mu=np.array([run.fermi]))

    s_doble = np.trace(run.sigma[0, 0]) / 3.0
    s_mitad = np.trace(mitad.sigma[0, 0]) / 3.0
    assert s_doble == pytest.approx(2.0 * s_mitad, rel=1e-12)
    # el número de Lorenz es un cociente: no puede moverse
    assert float(tr.lorenz(mitad, 0)[0]) == pytest.approx(
        float(tr.lorenz(run, 0)[0]), rel=1e-12)
