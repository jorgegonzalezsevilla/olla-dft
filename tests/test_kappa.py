"""Conductividad térmica de red.

Las pruebas que tocan phono3py se saltan si no está instalado: es una
dependencia opcional, y el resto del módulo (informe, acumulada, ajustes)
tiene que funcionar y probarse sin ella.
"""

import importlib.util

import numpy as np
import pytest

from qekit.core.errors import ErrorDeUso, FaltanDatos
from qekit.modules import kappa as K

# phono3py es un extra opcional y CI solo instala `.[test]`. Cuando este
# importorskip estaba a nivel de MÓDULO se saltaban las 30 pruebas, incluidas
# las del post-proceso, que no lo necesitan para nada: el módulo quedaba con
# 1 de 13 funciones ejercitadas y tres pruebas llevaban rotas desde la
# traducción al inglés sin que nadie lo viera. Ahora solo lo llevan las que
# tocan la librería.
necesita_phono3py = pytest.mark.skipif(
    importlib.util.find_spec("phono3py") is None,
    reason="phono3py es un extra opcional (pip install olla-dft[kappa])")


def _si():
    from ase.build import bulk
    return bulk("Si", "diamond", 5.43)


# ----------------------------------------------------------------------
# Configuraciones desplazadas
# ----------------------------------------------------------------------
@necesita_phono3py
def test_el_silicio_2x2x2_pide_57_configuraciones():
    """Es un número fijo: lo fija la simetría, no una elección.

    Si cambia, o cambió phono3py o cambió la detección de simetría, y en
    cualquiera de los dos casos hay que enterarse.
    """
    ph = K.preparar(_si(), (2, 2, 2))
    s3, s2 = K.configuraciones(ph)
    assert len(s3) == 57
    assert len(s3[0]) == 16
    assert s2 == []


@necesita_phono3py
def test_una_supercelda_aparte_para_la_parte_armonica():
    ph = K.preparar(_si(), (2, 2, 2), dim_fc2=(3, 3, 3))
    s3, s2 = K.configuraciones(ph)
    assert len(s3) == 57 and len(s3[0]) == 16
    assert len(s2) >= 1 and len(s2[0]) == 54


@necesita_phono3py
def test_las_configuraciones_son_estructuras_de_ase_periodicas():
    ph = K.preparar(_si(), (2, 2, 2))
    s3, _ = K.configuraciones(ph)
    a = s3[0]
    assert all(a.pbc)
    assert a.get_chemical_symbols().count("Si") == 16
    # y están DESPLAZADAS respecto de la supercelda perfecta
    b = s3[1]
    assert not np.allclose(a.get_positions(), b.get_positions())


@necesita_phono3py
def test_el_desplazamiento_es_el_pedido():
    """Se mide sobre las propias estructuras, no sobre la API de phono3py.

    Los nombres internos de phono3py cambian entre versiones; la distancia
    que separa la supercelda desplazada de la perfecta, no.
    """
    ph = K.preparar(_si(), (2, 2, 2), distancia=0.05)
    perfecta = K._ase(ph.supercell)
    d = (K.configuraciones(ph)[0][0].get_positions()
         - perfecta.get_positions())
    movidos = np.linalg.norm(d, axis=1)
    assert movidos.max() == pytest.approx(0.05, abs=1e-9)
    assert (movidos > 1e-9).sum() == 1        # solo uno se mueve en la fc3


@necesita_phono3py
def test_escribir_inputs_hace_una_carpeta_por_configuracion(tmp_path):
    from qekit.modules import sweep
    ph = K.preparar(_si(), (2, 2, 2))
    s3, _ = K.configuraciones(ph)
    s3 = s3[:4]
    common = sweep.prepare_common(s3[0], "/usr/share/espresso/pseudo",
                                  25, 100, True)
    carp = K.escribir_inputs(s3, tmp_path, common, kspacing=0.4)
    assert len(carp) == 4
    for d in carp:
        assert (d / "pw.in").exists()
        txt = (d / "pw.in").read_text()
        assert "calculation" in txt and "K_POINTS" in txt
        assert txt.count("Si ") >= 16 or txt.count("Si") >= 16


def test_leer_fuerzas_dice_cuales_faltan(tmp_path):
    d1, d2 = tmp_path / "d0000", tmp_path / "d0001"
    d1.mkdir(); d2.mkdir()
    with pytest.raises(FaltanDatos) as e:
        K.leer_fuerzas([d1, d2], 16)
    assert "d0000" in str(e.value) and "d0001" in str(e.value)


# ----------------------------------------------------------------------
# Post-proceso, sin necesidad de phono3py
# ----------------------------------------------------------------------
def _run_sintetico(n=1.0, nT=8):
    T = np.linspace(100, 800, nT)
    k = 100.0 * (T / 300.0) ** (-n)
    kap = np.zeros((nT, 6))
    kap[:, 0] = kap[:, 1] = kap[:, 2] = k
    return K.KappaRun(formula="Si2", dim=(3, 3, 3), temperaturas=T,
                      kappa=kap, i300=int(np.argmin(abs(T - 300))),
                      fuente="Quantum ESPRESSO", n_config=57, n_atomos=16)


def test_la_media_del_tensor_es_la_traza_entre_tres():
    run = _run_sintetico()
    run.kappa[:, 1] *= 2.0
    esperado = (run.kappa[:, 0] + run.kappa[:, 1] + run.kappa[:, 2]) / 3
    assert np.allclose(run.kappa_media, esperado)


@pytest.mark.parametrize("n", [0.5, 1.0, 1.5, 2.0])
def test_el_exponente_de_la_temperatura_se_recupera(n):
    run = _run_sintetico(n)
    assert K.exponente_temperatura(run) == pytest.approx(n, abs=1e-6)


def test_sin_datos_no_hay_exponente():
    run = K.KappaRun()
    assert K.exponente_temperatura(run) is None


def test_el_informe_dice_que_es_rta_y_que_faltan_los_isotopos():
    r = K.report(_run_sintetico())
    assert "RTA" in r
    assert "isotope" in r
    assert "four-phonon" in r      # los procesos de cuatro fonones


def test_el_informe_avisa_de_que_el_potencial_no_es_dft():
    run = _run_sintetico()
    run.fuente = "MACE"
    r = K.report(run)
    assert "not from DFT" in r
    assert "MACE" in r


def test_el_informe_avisa_de_una_supercelda_pequena():
    run = _run_sintetico()
    run.dim = (2, 2, 2)
    assert "supercell is 2×2×2, which is small" in K.report(run)


def test_una_supercelda_grande_no_dispara_el_aviso():
    """Antes buscaba «pequeña» en un informe en inglés: pasaba sin comprobar nada."""
    run = _run_sintetico()
    run.dim = (3, 3, 3)
    assert "which is small" not in K.report(run)


def test_el_informe_no_duplica_los_avisos_al_llamarlo_dos_veces():
    """Regresión: `report` hacía append sobre run.avisos.

    El CLI lo llama y después llama a `export`, que vuelve a llamarlo, así que
    cada WARNING salía dos veces en KAPPA.txt (y tres si algo más lo pedía).
    """
    run = _run_sintetico()
    run.fuente, run.dim = "MACE", (2, 2, 2)
    primero = K.report(run)
    assert primero.count("WARNING:") == 2
    assert K.report(run).count("WARNING:") == 2
    assert run.avisos == [], "report() no debe modificar el KappaRun"


def test_el_informe_reconoce_el_uno_partido_por_t_de_umklapp():
    assert "Umklapp" in K.report(_run_sintetico(1.0))


def test_export_escribe_kappa_y_el_informe(tmp_path):
    f = K.export(_run_sintetico(), str(tmp_path))
    assert any(x.endswith("KAPPA.dat") for x in f)
    assert any(x.endswith("KAPPA.txt") for x in f)
    dat = np.loadtxt(tmp_path / "KAPPA.dat")
    assert dat.shape == (8, 8)


# ----------------------------------------------------------------------
# Acumulada frente al recorrido libre medio
# ----------------------------------------------------------------------
def _run_con_modos(nq=40, nb=6, semilla=0):
    rng = np.random.default_rng(semilla)
    run = _run_sintetico()
    nT = len(run.temperaturas)
    run.pesos = np.ones(nq)
    run.velocidades = rng.random((nq, nb, 3)) * 50.0        # Å/ps
    run.gamma = rng.random((nT, nq, nb)) * 0.5 + 0.01       # THz
    run.cv = rng.random((nT, nq, nb)) + 0.1
    return run


def test_la_acumulada_va_de_cero_a_uno_y_no_baja():
    L, a = K.acumulada(_run_con_modos())
    assert a[-1] == pytest.approx(1.0)
    assert np.all(np.diff(a) >= -1e-12)
    assert np.all(np.diff(L) >= -1e-12)


def test_el_recorrido_al_noventa_por_ciento_es_mayor_que_al_cincuenta():
    run = _run_con_modos()
    assert (K.recorrido_representativo(run, 0.9)
            > K.recorrido_representativo(run, 0.5))


def test_sin_gammas_no_hay_acumulada():
    L, a = K.acumulada(_run_sintetico())
    assert L is None and a is None


def test_export_incluye_la_acumulada_cuando_la_hay(tmp_path):
    f = K.export(_run_con_modos(), str(tmp_path))
    assert any("recorrido" in x for x in f)
    d = np.loadtxt(tmp_path / "KAPPA_recorrido.dat")
    assert d.shape[1] == 2 and 0.0 <= d[:, 1].min() and d[:, 1].max() <= 1.0


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------
@pytest.mark.parametrize("txt,esp", [
    ("300", [300.0]),
    ("300,500", [300.0, 500.0]),
    ("100:300:3", [100.0, 200.0, 300.0]),
    ("300;500", [300.0, 500.0]),
])
def test_lista_de_temperaturas(txt, esp):
    from qekit.cli import _temperaturas
    assert _temperaturas(txt) == pytest.approx(esp)


@pytest.mark.parametrize("mal", ["100:800", "a:b:c", "100:800:0", "abc"])
def test_lista_de_temperaturas_mal_escrita(mal):
    from qekit.cli import _temperaturas
    with pytest.raises(ErrorDeUso):
        _temperaturas(mal)


# ----------------------------------------------------------------------
# Lo que se le pide a phono3py y con qué convenio se lee lo que devuelve.
# Dos fallos que no cambiaban nada visible: la opción quedaba inerte y el
# eje Λ salía escalado. Ninguna de las dos pruebas necesita la librería.
# ----------------------------------------------------------------------
class _Phono3pyFalso:
    """Registra lo que `resolver` le pide, sin calcular nada."""

    def __init__(self):
        self.kw = {}
        self.mesh_numbers = None
        self.forces = None
        self.phonon_forces = None
        self.thermal_conductivity = "tc"
        self.llamadas = []

    def __getattr__(self, nombre):
        if nombre.startswith(("produce_", "symmetrize_", "init_")):
            return lambda *a, **k: self.llamadas.append(nombre)
        raise AttributeError(nombre)

    def run_thermal_conductivity(self, **kw):
        self.kw = kw


def test_el_tamano_de_grano_se_pasa_en_micrometros():
    """Regresión: se convertía a Å (×1e4) y phono3py lo espera en µm.

    El propio phono3py calcula Γ_b = |v|·1e6·Å/(4π·boundary_mfp) y lo
    documenta en µm. Con el ×1e4 la dispersión por fronteras salía 10⁴ veces
    más débil: `--grain 0.1` daba la κ del cristal infinito mientras el
    informe decía que se había aplicado el grano.
    """
    ph = _Phono3pyFalso()
    K.resolver(ph, np.zeros((1, 2, 3)), frontera_um=0.1, malla=5)
    assert ph.kw["boundary_mfp"] == pytest.approx(0.1)


def test_sin_grano_se_pide_el_cristal_infinito():
    ph = _Phono3pyFalso()
    K.resolver(ph, np.zeros((1, 2, 3)), malla=5)
    assert ph.kw["boundary_mfp"] >= 1e6      # 1 m: sin fronteras
    assert ph.mesh_numbers == [5, 5, 5]


def _run_de_un_modo(gamma=0.5, v=200.0):
    """Un único modo, para poder comprobar Λ a mano."""
    run = K.KappaRun()
    run.temperaturas = np.array([300.0])
    run.i300 = 0
    run.gamma = np.array([[[gamma]]])                  # (nT, nq, nb)
    run.velocidades = np.array([[[v, 0.0, 0.0]]])      # (nq, nb, 3)
    run.cv = np.array([[[1.0]]])                       # (nT, nq, nb)
    run.pesos = np.array([1.0])
    return run


def test_la_vida_media_lleva_el_2pi_del_convenio_de_phono3py():
    """Λ = |v|/(2·2π·Γ), no |v|/(2Γ).

    Un 2 pasa Γ de HWHM a anchura total y el 2π pasa de frecuencia cíclica a
    angular, porque la Γ que devuelve phono3py viene en THz ordinarios. Sin el
    2π todos los recorridos libres medios salían 6.28 veces más largos.
    """
    gamma, v = 0.5, 200.0
    L, _ = K.acumulada(_run_de_un_modo(gamma, v))
    assert L[0] == pytest.approx(v / (2.0 * 2.0 * np.pi * gamma))


@necesita_phono3py
def test_el_recorrido_coincide_con_el_de_phono3py():
    """La comprobación que zanja el convenio: contra la propia librería."""
    from phono3py.conductivity.utils import get_mfp

    gamma, v = 0.5, 200.0
    L, _ = K.acumulada(_run_de_un_modo(gamma, v))
    suyo = get_mfp(np.array([[gamma]]), np.array([[[v, 0.0, 0.0]]]))
    assert L[0] == pytest.approx(float(suyo[0, 0]))


def test_la_fraccion_acumulada_no_depende_del_convenio():
    """El 2π se cancela al normalizar: solo se mueve el eje Λ, no la curva."""
    run = K.KappaRun()
    run.temperaturas, run.i300 = np.array([300.0]), 0
    run.gamma = np.array([[[0.5, 2.0]]])
    run.velocidades = np.array([[[100.0, 0, 0], [300.0, 0, 0]]])
    run.cv = np.array([[[1.0, 1.0]]])
    run.pesos = np.array([1.0])
    L, a = K.acumulada(run)
    assert a[-1] == pytest.approx(1.0)
    assert np.all(np.diff(a) >= -1e-12)


def test_recoger_pasa_lo_de_phono3py_a_los_campos_del_run():
    """Sin la librería: basta con un objeto que exponga los mismos atributos."""
    from types import SimpleNamespace

    T = np.array([100.0, 300.0, 500.0])
    tc = SimpleNamespace(
        temperatures=T,
        kappa=np.arange(3 * 6, dtype=float).reshape(1, 3, 6),
        frequencies=np.ones((4, 2)),
        grid_weights=np.ones(4),
        gamma=np.full((1, 3, 4, 2), 0.5),
        group_velocities=np.full((4, 2, 3), 100.0),
        mode_heat_capacities=np.ones((3, 4, 2)),
    )
    run = K.recoger(K.KappaRun(), None, tc, (5, 5, 5))
    assert run.malla == (5, 5, 5)
    assert run.i300 == 1, "i300 tiene que ser el índice del T más cercano a 300"
    assert run.kappa.shape == (3, 6)
    assert run.gamma.shape == (3, 4, 2)
    assert run.avisos == []


def test_recoger_avisa_si_faltan_los_datos_por_modo():
    """Regresión: un `except Exception: pass` los borraba sin decir nada.

    κ(T) sigue saliendo, pero la sección de recorrido libre medio desaparecía
    junto con su fichero y su figura, y el usuario no sabía por qué.
    """
    from types import SimpleNamespace

    tc = SimpleNamespace(
        temperatures=np.array([300.0]),
        kappa=np.zeros((1, 1, 6)),
        frequencies=np.ones((2, 1)),
        grid_weights=np.ones(2),
    )                              # sin gamma / group_velocities / cv
    run = K.recoger(K.KappaRun(), None, tc, (5, 5, 5))
    assert run.kappa is not None
    assert run.gamma is None
    assert any("mean free path" in a for a in run.avisos)


def test_plot_escribe_las_dos_figuras_cuando_hay_acumulada(tmp_path):
    escritos = K.plot(_run_con_modos(), str(tmp_path / "k"), formats="png")
    assert escritos
    assert any("recorrido" in f for f in escritos), \
        f"falta la figura del recorrido libre medio: {escritos}"


def test_plot_sin_kappa_lo_dice(tmp_path):
    with pytest.raises(FaltanDatos):
        K.plot(K.KappaRun(), str(tmp_path / "k"), formats="png")


# ----------------------------------------------------------------------
# La acumulada tiene que cerrar con la κ que se informa arriba
# ----------------------------------------------------------------------
def _run_dos_modos(mode_kappa=True, iso=None, frontera=None):
    """Dos modos con κ conocida por modo, para comprobar el reparto."""
    run = K.KappaRun()
    run.temperaturas, run.i300 = np.array([300.0]), 0
    run.gamma = np.array([[[1.0, 4.0]]])                   # (nT, nq, nb)
    run.velocidades = np.array([[[300.0, 0, 0], [100.0, 0, 0]]])
    run.cv = np.array([[[1.0, 1.0]]])
    run.pesos = np.array([1.0])
    run.gamma_iso = iso
    run.frontera = frontera
    if mode_kappa:
        # (nT, nq, nb, 6) en Voigt: 30 y 10 W/mK de traza/3
        run.mode_kappa = np.zeros((1, 1, 2, 6))
        run.mode_kappa[0, 0, 0, :3] = 30.0
        run.mode_kappa[0, 0, 1, :3] = 10.0
    return run


def test_la_acumulada_reparte_la_kappa_de_cada_modo():
    """El peso es `mode_kappa`, la descomposición que phono3py ya trae.

    Con 30 y 10 W/mK, el modo pequeño (Λ menor) tiene que llevar el 25 % y el
    grande el 75 % restante. Reconstruir el peso a mano como C·v²·τ daba otra
    curva en cuanto había más de un canal de dispersión.
    """
    L, a = K.acumulada(_run_dos_modos())
    assert len(a) == 2
    assert L[0] < L[1], "el modo lento y ancho tiene que ir primero"
    assert a[0] == pytest.approx(0.25)
    assert a[1] == pytest.approx(1.0)


def test_sin_mode_kappa_se_recae_en_la_reconstruccion():
    """Versiones viejas de phono3py no la exponen: la curva sigue saliendo."""
    L, a = K.acumulada(_run_dos_modos(mode_kappa=False))
    assert L is not None and a[-1] == pytest.approx(1.0)


def test_gamma_total_suma_isotopos_y_fronteras():
    """Regresión: Λ usaba solo la parte fonón-fonón.

    Con `--isotopes` o `--grain`, phono3py suma esos canales para calcular κ,
    así que el informe daba la κ del cristal con isótopos y grano y, debajo, un
    recorrido libre medio del cristal puro e infinito, sin decirlo.
    """
    solo_phph = K.gamma_total(_run_dos_modos(), 0)
    assert np.allclose(solo_phph, [[1.0, 4.0]])

    con_iso = K.gamma_total(_run_dos_modos(iso=np.array([[0.5, 0.5]])), 0)
    assert np.allclose(con_iso, [[1.5, 4.5]])

    # Γ_b = |v|·1e6·Å/(4π·L), la fórmula literal de phono3py
    run = _run_dos_modos(frontera=2.0)
    esperado = np.array([[1.0, 4.0]]) + np.array([300.0, 100.0]) * 1e6 * 1e-10 / (
        4.0 * np.pi * 2.0)
    assert np.allclose(K.gamma_total(run, 0), esperado)


def test_un_grano_acorta_todos_los_recorridos():
    """Ningún fonón puede recorrer más que el grano, y antes sí lo hacía."""
    sin_grano = K.recorrido_representativo(_run_dos_modos(), 0.9)
    con_grano = K.recorrido_representativo(_run_dos_modos(frontera=0.05), 0.9)
    assert con_grano < sin_grano
    assert con_grano < 0.05 * 1e4, "Λ90 tiene que quedar por debajo del grano"


def _run_para_informe(malla):
    run = K.KappaRun()
    run.temperaturas, run.i300 = np.array([300.0]), 0
    run.kappa = np.full((1, 6), 100.0)
    run.dim, run.formula, run.malla = (3, 3, 3), "Si2", malla
    run.fuente = "Quantum ESPRESSO"
    run.gamma = np.array([[[1.0, 4.0]]])
    run.velocidades = np.array([[[300.0, 0, 0], [100.0, 0, 0]]])
    run.cv = np.array([[[1.0, 1.0]]])
    run.pesos = np.array([1.0])
    run.mode_kappa = np.zeros((1, 1, 2, 6))
    run.mode_kappa[0, 0, 0, :3] = 30.0
    run.mode_kappa[0, 0, 1, :3] = 10.0
    return run


def test_una_malla_floja_no_sirve_para_el_noventa_por_ciento():
    """Λ90 vive en la cola de recorridos largos, que una malla floja no muestrea.

    Medido en Si con Stillinger-Weber (supercelda 3×3×3): de 11³ a 31³, Λ50 va
    de 0.68 a 0.82 µm (+21 %) mientras Λ90 pasa de 4.4 a 21.9 µm, cinco veces.
    Citar Λ90 con la malla por omisión, o dimensionar un grano con él, no vale.
    """
    assert "NOT for Λ90" in K.report(_run_para_informe((13, 13, 13)))
    assert "NOT for Λ90" not in K.report(_run_para_informe((31, 31, 31)))
