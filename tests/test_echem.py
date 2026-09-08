"""Electrodo de hidrógeno computacional."""

import numpy as np
import pytest

from qekit.core.errors import ErrorDeUso
from qekit.modules import echem


# ----------------------------------------------------------------------
# HER
# ----------------------------------------------------------------------
def test_el_platino_esta_en_la_cumbre_del_volcan():
    """Pt(111): E_ads(H) = -0.33 eV y ZPE-TdS = +0.24 dan dG_H = -0.09."""
    e = echem.her(-0.33)
    assert e.dG_H == pytest.approx(-0.09, abs=1e-9)
    assert 'top of the volcano' in echem.report(e)


def test_los_dos_pasos_de_la_her_son_opuestos():
    e = echem.her(-0.33)
    assert e.pasos[0][1] == pytest.approx(-e.pasos[1][1])


def test_el_sobrepotencial_de_la_her_es_el_paso_peor():
    e = echem.her(-0.33)
    assert e.sobrepotencial == pytest.approx(0.09, abs=1e-9)
    e2 = echem.her(0.06)          # dG_H = +0.30
    assert e2.sobrepotencial == pytest.approx(0.30, abs=1e-9)


def test_se_reconoce_cada_rama_del_volcan():
    assert 'binds too strongly' in echem.report(echem.her(-1.0))
    assert 'barely adsorbs' in echem.report(echem.her(+1.0))


def test_la_correccion_se_puede_dar_a_mano():
    assert echem.her(-0.33, correccion=0.0).dG_H == pytest.approx(-0.33)


def test_avisa_de_que_la_correccion_es_de_tabla():
    assert 'not one computed' in echem.report(echem.her(-0.33))
    assert 'not one computed' not in echem.report(echem.her(-0.33, correccion=0.2))


# ----------------------------------------------------------------------
# OER
# ----------------------------------------------------------------------
SIN_CORR = {"OH": 0, "O": 0, "OOH": 0}
RUO2 = {"OH": 0.77, "O": 2.16, "OOH": 3.87}


def test_el_ruo2_da_su_sobrepotencial_publicado():
    """Man et al. 2011: RuO2(110) tiene eta ~ 0.48 V."""
    e = echem.oer(RUO2, correcciones=SIN_CORR)
    assert e.sobrepotencial == pytest.approx(0.481, abs=0.01)


def test_los_cuatro_pasos_suman_el_total_del_agua():
    """Por construccion: el cuarto sale por diferencia con 4.92 eV."""
    e = echem.oer(RUO2, correcciones=SIN_CORR)
    assert sum(g for _, g in e.pasos) == pytest.approx(echem.DG_AGUA_TOTAL,
                                                       abs=1e-9)


def test_el_paso_limitante_del_ruo2_es_la_formacion_de_OOH():
    e = echem.oer(RUO2, correcciones=SIN_CORR)
    nombre, _ = e.limitante
    assert "OOH" in nombre


def test_faltan_intermedios():
    with pytest.raises(ErrorDeUso) as exc:
        echem.oer({"OH": 0.8})
    assert "O" in str(exc.value) and "OOH" in str(exc.value)


def test_avisa_si_el_cuarto_paso_sale_negativo():
    """Si los tres calculados ya suman mas de 4.92, algo falla."""
    e = echem.oer({"OH": 1.5, "O": 3.0, "OOH": 5.0}, correcciones=SIN_CORR)
    assert any('NEGATIVE' in a for a in e.avisos)


def test_la_relacion_de_escala_se_comprueba():
    e = echem.oer(RUO2, correcciones=SIN_CORR)
    txt = echem.report(e)
    assert "3.100" in txt and "3.2" in txt


def test_una_escala_rara_se_senala():
    e = echem.oer({"OH": 0.5, "O": 2.0, "OOH": 5.0}, correcciones=SIN_CORR)
    assert 'before celebrating' in echem.report(e)


# ----------------------------------------------------------------------
# potencial y pH
# ----------------------------------------------------------------------
def test_el_potencial_baja_todos_los_pasos_por_igual():
    e = echem.oer(RUO2, correcciones=SIN_CORR)
    g0 = [g for _, g in e.dG(0.0, 0.0)]
    g1 = [g for _, g in e.dG(1.0, 0.0)]
    assert all(b == pytest.approx(a - 1.0) for a, b in zip(g0, g1))


def test_al_potencial_limitante_ningun_paso_es_cuesta_arriba():
    e = echem.oer(RUO2, correcciones=SIN_CORR)
    gs = [g for _, g in e.dG(e.U_limitante, 0.0)]
    assert max(gs) == pytest.approx(0.0, abs=1e-9)


def test_el_ph_entra_con_kT_ln10():
    """El tamaño del término de pH es kT·ln10 por unidad, y el SIGNO depende
    de la reacción: la OER suelta H⁺ (subir el pH la ayuda) y la HER lo
    consume (subir el pH la estorba)."""
    magnitud = 7.0 * echem.KB_EV * 298.15 * np.log(10.0)

    o = echem.oer(RUO2, correcciones=SIN_CORR)
    assert o.dG(0.0, 0.0)[0][1] - o.dG(0.0, 7.0)[0][1] == pytest.approx(
        magnitud, rel=1e-9)

    h = echem.her(-0.33)
    assert h.dG(0.0, 7.0)[0][1] - h.dG(0.0, 0.0)[0][1] == pytest.approx(
        magnitud, rel=1e-9), \
        "en la HER, subir el pH tiene que ENCARECER el paso, no abaratarlo"


def test_la_her_necesita_potencial_catodico():
    """Regresión: la HER llevaba el signo de potencial de la OER.

    Sus pasos consumen H⁺+e⁻, así que ΔG = ΔG₀ + eU: el potencial que pone
    a cero el paso peor es NEGATIVO frente al RHE. Con el signo de la OER se
    informaba U_L = +|ΔG| para una reacción catódica, y la escalera de
    energías se movía al revés al aplicar potencial.
    """
    e = echem.her(0.0)                       # ΔG_H* = +0.24 eV
    assert e.dG_H == pytest.approx(0.24)
    assert e.U_limitante == pytest.approx(-0.24), \
        "la HER es catódica: su potencial limitante es negativo"
    # η se sigue informando positivo (es lo que hay que aplicar de más)
    assert e.sobrepotencial == pytest.approx(0.24)
    # y a ese potencial ningún paso queda cuesta arriba
    assert max(g for _, g in e.dG(e.U_limitante, 0.0)) == pytest.approx(
        0.0, abs=1e-9)


def test_el_pt_sigue_en_la_cumbre_del_volcan():
    """El descriptor no cambia con el signo: |ΔG_H*| de Pt(111) es 0.09 eV."""
    e = echem.her(-0.33)
    assert e.dG_H == pytest.approx(-0.09)
    assert e.sobrepotencial == pytest.approx(0.09)
    assert e.U_limitante == pytest.approx(-0.09)


def test_el_diagrama_de_pourbaix_tiene_la_forma_pedida():
    e = echem.oer(RUO2, correcciones=SIN_CORR)
    d = echem.pourbaix(e, U=np.linspace(0, 2, 11), pH=np.linspace(0, 14, 8))
    assert d["dG_limitante"].shape == (8, 11)
    # sube con el pH y baja con el potencial
    assert d["dG_limitante"][-1, 0] < d["dG_limitante"][0, 0]
    assert d["dG_limitante"][0, -1] < d["dG_limitante"][0, 0]


def test_export(tmp_path):
    e = echem.her(-0.33)
    f = echem.export(e, str(tmp_path))
    assert len(f) == 2
    assert "sobrepotencial_V" in (tmp_path / "ECHEM.dat").read_text()
