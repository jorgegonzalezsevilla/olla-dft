"""Los tres módulos que la suite no llegaba a ejecutar nunca.

`health` (lo que corre `olla-dft doctor --system`), `converge` y `exfoliate`
tenían cero funciones ejercitadas: se importaban y como mucho se instanciaba
su dataclass. Eran 628 líneas sin una sola aserción encima.
"""

import json

import numpy as np
import pytest
from ase import Atoms
from ase.build import bulk

from qekit.core.errors import ErrorDeUso, FaltanDatos
from qekit.modules import converge, exfoliate, health

SIN_PSEUDOS = "/no/existe/pseudos"


def grafito():
    """Dos capas de grafeno apiladas: el caso laminar de libro."""
    a, c = 2.46, 6.70
    cell = [[a, 0.0, 0.0], [-a / 2, a * np.sqrt(3) / 2, 0.0], [0.0, 0.0, c]]
    return Atoms("C4", cell=cell, pbc=True, scaled_positions=[
        (0.0, 0.0, 0.0), (1 / 3, 2 / 3, 0.0),
        (0.0, 0.0, 0.5), (2 / 3, 1 / 3, 0.5)])


# ----------------------------------------------------------------------
# health: el diagnóstico de la instalación
# ----------------------------------------------------------------------
def test_el_diagnostico_encuentra_el_paquete_y_sus_dependencias():
    r = health.check()
    por_codigo = {c["code"]: c for c in r["checks"]}
    assert por_codigo["qekit.version"]["level"] == "ok"
    assert por_codigo["python.version"]["level"] == "ok"
    dep = por_codigo["python.dependencies"]
    assert dep["level"] == "ok", f"faltan dependencias: {dep['detail']}"
    # la evidencia es un JSON con la versión de cada una
    versiones = json.loads(dep["evidence"])
    assert set(versiones) == set(health.DEPENDENCIES)


def test_las_dependencias_declaradas_son_las_del_pyproject():
    """Si pyproject gana una dependencia y health no, el diagnóstico miente."""
    import importlib.metadata
    import re
    requeridas = set()
    for linea in importlib.metadata.requires("olla-dft") or []:
        if "extra ==" in linea:          # los extras son opcionales
            continue
        m = re.match(r"^\s*([A-Za-z0-9_.\-]+)", linea)
        if m:
            requeridas.add(m.group(1).lower())
    assert requeridas == {d.lower() for d in health.DEPENDENCIES}


def test_el_diagnostico_no_da_por_buena_una_instalacion_sin_pw(monkeypatch):
    """Regresión: tener ph.x pero no pw.x se informaba como 'ok'.

    Sin pw.x no se puede correr NINGÚN cálculo, así que decir que Quantum
    ESPRESSO está bien porque aparece otro binario engaña.
    """
    monkeypatch.setattr(health.shutil, "which",
                        lambda n: "/usr/bin/ph.x" if n.startswith("ph") else None)
    qe = {c["code"]: c for c in health.check()["checks"]}["qe.binaries"]
    assert qe["level"] == "warn"
    assert "NOT pw.x" in qe["detail"]


def test_el_diagnostico_busca_el_nombre_del_binario_de_cada_sistema(monkeypatch):
    """En Windows los binarios son pw.exe; buscar solo pw.x no encuentra nada."""
    vistos = []

    def which(nombre):
        vistos.append(nombre)
        return None

    monkeypatch.setattr(health.shutil, "which", which)
    health.check()
    assert any(n.endswith(".exe") for n in vistos), \
        "solo se probaron nombres .x: en Windows no encontraría QE"


def test_la_memoria_se_mide_en_los_tres_sistemas(monkeypatch):
    """Regresión: se leía solo /proc/meminfo, así que fuera de Linux no había medida."""
    from qekit.core import runner
    monkeypatch.setattr(runner, "memoria_libre_gb", lambda: 7.5)
    mem = {c["code"]: c for c in health.check()["checks"]}["resources.memory"]
    assert mem["level"] == "ok"
    assert "7.50" in mem["detail"]


def test_el_diagnostico_avisa_de_la_memoria_escasa(monkeypatch):
    from qekit.core import runner
    monkeypatch.setattr(runner, "memoria_libre_gb", lambda: 0.2)
    mem = {c["code"]: c for c in health.check()["checks"]}["resources.memory"]
    assert mem["level"] == "fail"


def test_un_fallo_bloquea_el_veredicto(monkeypatch):
    monkeypatch.setattr(health.sys, "version_info", (3, 8, 0))
    r = health.check()
    assert r["ok"] is False and r["fails"] >= 1
    assert "BLOCKED" in health.report(r)


def test_el_informe_lista_todas_las_comprobaciones():
    r = health.check()
    texto = health.report(r)
    assert "READY" in texto or "BLOCKED" in texto
    for c in r["checks"]:
        assert c["title"] in texto


def test_un_proyecto_que_no_existe_se_informa_sin_reventar(tmp_path):
    r = health.check(project_path=str(tmp_path / "no_hay"))
    fallo = {c["code"]: c for c in r["checks"]}["project.load"]
    assert fallo["level"] == "fail" and r["ok"] is False


# ----------------------------------------------------------------------
# converge: análisis de la serie
# ----------------------------------------------------------------------
def _serie(energias, natoms=2, umbral=1.0, kind="ecutwfc"):
    r = converge.ConvergenceRun(kind=kind, natoms=natoms, threshold=umbral)
    r.values = [30.0 + 10 * i for i in range(len(energias))]
    r.labels = [f"ecutwfc = {v:g} Ry" for v in r.values]
    r.energies = list(energias)
    return r


def test_la_diferencia_se_mide_contra_el_punto_mas_denso():
    """El criterio del módulo: referencia = el último, no el anterior."""
    r = _serie([-1.0, -1.002, -1.0030, -1.0031])
    d = r.per_atom_diffs()
    # (E - E_ref) en meV/átomo, con 2 átomos
    assert d[-1] == pytest.approx(0.0)
    assert d[0] == pytest.approx(abs(-1.0 + 1.0031) * 1000 / 2)


def test_converge_donde_la_cola_entera_baja_del_umbral():
    # 0.05 y 0.0 meV/átomo al final; los dos primeros muy lejos
    r = _serie([-1.0, -1.5, -1.9999, -2.0], umbral=1.0)
    assert r.converged_index() == 2
    assert "CONVERGES at" in converge.report(r)


def test_no_converge_si_el_ultimo_salto_es_grande():
    r = _serie([-1.0, -1.5, -2.0, -3.0], umbral=1.0)
    texto = converge.report(r)
    assert "NOT converged" in texto
    assert "Extend the series" in texto


def test_el_ultimo_punto_solo_no_demuestra_convergencia():
    """Regresión: el último es su propia referencia y cumple siempre.

    Por eso `converged_index` no puede devolver None con datos suficientes, y
    la rama "NOT converged" que había aparte era inalcanzable: el informe
    prometía un mensaje que nadie llegaba a ver nunca.
    """
    r = _serie([-1.0, -2.0, -3.0], umbral=1.0)
    assert r.converged_index() == len(r.values) - 1
    assert "no margin" in converge.report(r)
    assert "NOT converged" in converge.report(r)


def test_un_punto_fallido_no_rompe_el_analisis():
    r = _serie([-1.0, None, -1.0001, -1.00011])
    d = r.per_atom_diffs()
    assert np.isnan(d[1])
    assert r.converged_index() is not None
    assert "FAILED" in converge.report(r)


def test_sin_puntos_suficientes_lo_dice_en_vez_de_calcular():
    r = _serie([-1.0, None, None])
    assert r.per_atom_diffs().size == 0
    assert r.converged_index() is None
    assert "Not enough finished calculations" in converge.report(r)


def test_el_informe_recomienda_el_valor_que_convergio():
    r = _serie([-1.0, -1.5, -1.9999, -2.0], umbral=1.0)
    assert "--ecutwfc 50" in converge.report(r)


def test_export_escribe_dat_y_txt_en_utf8(tmp_path):
    r = _serie([-1.0, -1.5, -1.9999, -2.0])
    escritos = converge.export(r, str(tmp_path))
    assert len(escritos) == 2
    dat = (tmp_path / "CONVERGENCIA.dat").read_text(encoding="utf-8")
    assert "meV/atom" in dat
    # el informe lleva ΔE: si se escribiera con la codificación de la locale,
    # en Windows esto reventaría al escribir, no al leer
    txt = (tmp_path / "CONVERGENCIA.txt").read_text(encoding="utf-8")
    assert "Δ" in txt


def test_la_curva_necesita_dos_puntos_terminados(tmp_path):
    """Regresión: con uno solo, matplotlib moría por longitudes distintas."""
    r = _serie([-1.0, None, None])
    with pytest.raises(FaltanDatos, match="TWO finished"):
        converge.plot(r, str(tmp_path / "c"), formats="png")


def test_la_curva_se_dibuja_con_puntos_suficientes(tmp_path):
    r = _serie([-1.0, -1.5, -1.9999, -2.0])
    escritos = converge.plot(r, str(tmp_path / "c"), formats="png")
    assert escritos and (tmp_path / "c.png").exists()


# ----------------------------------------------------------------------
# converge: preparación de la serie
# ----------------------------------------------------------------------
def test_prepara_la_serie_de_ecutwfc(tmp_path):
    r, rep = converge.prepare(bulk("Si", "diamond", 5.43), "ecutwfc",
                              outdir=str(tmp_path), pseudo_dir=SIN_PSEUDOS)
    assert len(r.jobs) == len(r.values) == 8
    assert r.values[0] == 30.0 and r.values[-1] == 100.0
    assert all(j.input_path.exists() for j in r.jobs)
    assert "pseudopotentials are missing" in rep


def test_prepara_la_serie_de_mallas_sin_repetir(tmp_path):
    r, _ = converge.prepare(bulk("Si", "diamond", 5.43), "kmesh",
                            outdir=str(tmp_path), pseudo_dir=SIN_PSEUDOS)
    mallas = [j.meta["grid"] for j in r.jobs]
    assert len(mallas) == len(set(mallas)), "hay mallas repetidas en la serie"
    assert r.values == sorted(r.values), "la serie no va de menos a más densa"


def test_un_tipo_de_convergencia_inventado_es_error_de_uso(tmp_path):
    with pytest.raises(ErrorDeUso, match="unknown convergence type"):
        converge.prepare(bulk("Si", "diamond", 5.43), "ecutfoo",
                         outdir=str(tmp_path), pseudo_dir=SIN_PSEUDOS)


@pytest.mark.parametrize("valor", ["4x4", "x", "8x8xz", "0x0x0", "hola", "-1"])
def test_una_malla_mal_escrita_es_error_de_uso(valor):
    """Regresión: '4x4' llegaba como tupla de dos y salía un IndexError."""
    with pytest.raises(ErrorDeUso):
        converge._grid_pedida(valor, bulk("Si", "diamond", 5.43))


def test_una_malla_bien_escrita_se_acepta_de_las_dos_formas():
    si = bulk("Si", "diamond", 5.43)
    assert converge._grid_pedida("8x8x8", si) == (8, 8, 8)
    assert converge._grid_pedida("0.20", si) == converge._grid_pedida(0.20, si)


def test_collect_recoge_las_energias_de_los_resultados(tmp_path):
    r, _ = converge.prepare(bulk("Si", "diamond", 5.43), "ecutwfc",
                            outdir=str(tmp_path), values=[30, 40],
                            pseudo_dir=SIN_PSEUDOS)

    class _Res:
        def __init__(self, job, ok, energy):
            self.job, self.ok, self.energy = job, ok, energy

    converge.collect(r, [_Res(r.jobs[0], True, -1.0),
                         _Res(r.jobs[1], False, None)])
    assert r.energies == [-1.0, None]


def test_collect_sin_resultados_lee_las_carpetas_y_no_revienta(tmp_path):
    r, _ = converge.prepare(bulk("Si", "diamond", 5.43), "ecutwfc",
                            outdir=str(tmp_path), values=[30, 40],
                            pseudo_dir=SIN_PSEUDOS)
    converge.collect(r)          # las carpetas están vacías: todo None
    assert r.energies == [None, None]


# ----------------------------------------------------------------------
# exfoliate
# ----------------------------------------------------------------------
def test_prepara_los_dos_calculos_del_grafito(tmp_path):
    r, rep = exfoliate.prepare(grafito(), outdir=str(tmp_path),
                               pseudo_dir=SIN_PSEUDOS)
    assert r.n_layers == 2
    assert r.natoms_bulk == 4 and r.natoms_slab == 2
    assert [j.meta["which"] for j in r.jobs] == ["bulk", "slab"]
    assert all(j.input_path.exists() for j in r.jobs)
    # área en el plano de la celda hexagonal: |a x b|
    assert r.area == pytest.approx(2.46 ** 2 * np.sqrt(3) / 2, rel=1e-9)


def test_la_monocapa_lleva_un_solo_punto_k_en_la_normal(tmp_path):
    exfoliate.prepare(grafito(), outdir=str(tmp_path),
                               pseudo_dir=SIN_PSEUDOS)
    entrada = (tmp_path / "monocapa" / "pw.in").read_text(encoding="utf-8")
    malla = entrada.split("K_POINTS")[1].split("\n")[1].split()
    assert malla[2] == "1", f"la normal no lleva 1 punto: {malla}"
    # y la malla en el plano se conserva igual que la del bulk
    assert malla[0] == malla[1] and int(malla[0]) > 1


def test_sin_correccion_de_dispersion_el_reporte_lo_advierte(tmp_path):
    _, rep = exfoliate.prepare(grafito(), outdir=str(tmp_path),
                               pseudo_dir=SIN_PSEUDOS)
    assert "NO van der Waals correction" in rep


def test_una_estructura_no_laminar_es_error_de_uso(tmp_path):
    with pytest.raises(ErrorDeUso, match="no layers"):
        exfoliate.prepare(bulk("Si", "diamond", 5.43), outdir=str(tmp_path),
                          pseudo_dir=SIN_PSEUDOS)


def test_la_energia_de_exfoliacion_sale_de_la_formula_del_modulo():
    r = exfoliate.ExfoliationRun(n_layers=2, area=5.24, natoms_bulk=4,
                                 natoms_slab=2, layer_formulas=["C2", "C2"])
    r.E_bulk, r.E_slab = -20.0, -9.9
    texto = exfoliate.report_result(r)
    esperado = (r.E_slab - r.E_bulk / 2) / r.area * exfoliate.EV_A2_TO_J_M2
    assert f"{esperado:8.4f} J/m²" in texto


def test_una_energia_negativa_se_advierte():
    r = exfoliate.ExfoliationRun(n_layers=2, area=5.24, natoms_bulk=4,
                                 natoms_slab=2, layer_formulas=["C2", "C2"])
    r.E_bulk, r.E_slab = -20.0, -10.5      # monocapa por debajo del bulk
    assert "NEGATIVE" in exfoliate.report_result(r)


def test_si_falta_un_calculo_lo_dice_en_vez_de_dividir():
    r = exfoliate.ExfoliationRun(n_layers=2, area=5.24, natoms_slab=2)
    r.E_bulk, r.E_slab = -20.0, None
    texto = exfoliate.report_result(r)
    assert "missing" in texto and "monolayer" in texto


def test_capas_distintas_invalidan_el_reparto_por_numero_de_capas(tmp_path):
    """Regresión: E_bulk/N solo vale si las capas son iguales.

    En una heteroestructura, `make_slab` aísla UNA capa y dividir entre N
    mezcla capas distintas. Antes el número salía igual, sin decir nada.
    """
    a, c = 2.46, 14.0
    cell = [[a, 0.0, 0.0], [-a / 2, a * np.sqrt(3) / 2, 0.0], [0.0, 0.0, c]]
    hetero = Atoms("CCBNN", cell=cell, pbc=True, positions=[
        (0, 0, 0), (1.2, 0.7, 0),
        (0, 0, 7.0), (1.2, 0.7, 7.0), (2.4, 1.4, 7.0)])

    r, rep = exfoliate.prepare(hetero, outdir=str(tmp_path),
                               pseudo_dir=SIN_PSEUDOS)
    assert r.equivalent_layers is False
    assert "NOT equivalent" in rep

    r.E_bulk, r.E_slab = -20.0, -9.9
    assert "NOT equivalent" in exfoliate.report_result(r)


def test_capas_iguales_no_disparan_el_aviso(tmp_path):
    r, rep = exfoliate.prepare(grafito(), outdir=str(tmp_path),
                               pseudo_dir=SIN_PSEUDOS)
    assert r.equivalent_layers is True
    assert "NOT equivalent" not in rep


def test_collect_reparte_las_energias_por_su_etiqueta(tmp_path):
    r, _ = exfoliate.prepare(grafito(), outdir=str(tmp_path),
                             pseudo_dir=SIN_PSEUDOS)

    class _Res:
        def __init__(self, job, ok, energy):
            self.job, self.ok, self.energy = job, ok, energy

    exfoliate.collect(r, [_Res(r.jobs[1], True, -9.9),
                          _Res(r.jobs[0], True, -20.0)])
    assert r.E_bulk == -20.0 and r.E_slab == -9.9
