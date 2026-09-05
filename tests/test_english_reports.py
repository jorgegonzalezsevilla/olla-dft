"""Human scientific reports are English; stored scientific contracts stay stable."""
import copy
from pathlib import Path

import numpy as np
import pytest

from qekit import cli
from qekit.core import i18n
from qekit.modules import datasheet, effmass, eos, project, quality, report


@pytest.mark.parametrize('language', ['en', 'es', 'de'])
def test_scientific_cli_output_is_english_in_every_interface(language, capsys):
    source = Path(__file__).parent / 'datos/Si_relajado.cif'
    assert cli.main(['--language', language, 'info', str(source)]) == 0
    output = capsys.readouterr().out
    assert 'Volume' in output and 'Å' in output
    assert 'Volumen' not in output


def test_eos_report_and_fitted_numbers_do_not_depend_on_language(monkeypatch):
    volumes = np.linspace(36., 44., 11)
    run = eos.EOSRun(volumes=volumes.tolist(), energies=eos.birch_murnaghan(
        volumes, -10., 40., .6, 4.2).tolist())
    run.fits = eos.fit_all(run)
    before = copy.deepcopy(run)
    texts = []
    for language in ('en', 'es', 'de'):
        monkeypatch.setattr(i18n, 'get_language', lambda: language)
        texts.append(eos.report(run))
    assert texts[0] == texts[1] == texts[2]
    assert 'Equation of state' in texts[0]
    assert run == before
    assert abs(run.fits[eos.DEFAULT_EQ].V0 - 40.) < 1e-5


def test_legacy_datasheet_keys_render_in_english_without_mutation():
    sheet = datasheet.Ficha(formula='Si', resultados={'Elásticas': [
        {'magnitud': 'C₁₁', 'valor': 159.912345, 'unidad': 'GPa', 'nota': ''}]},
        parametros={'funcional': 'PBE', 'malla_k': '8x8x8', 'ocupaciones': 'fixed'})
    original = copy.deepcopy(sheet)
    for content in (datasheet.markdown(sheet), datasheet.html(sheet)):
        assert 'Elastic' in content and 'functional' in content
        assert 'Methods' in content and '159.912345' in content
        assert 'Elásticas' not in content
    assert sheet == original


def test_effective_mass_report_preserves_legacy_carriers_and_directions(tmp_path):
    run = effmass.EffMassRun(fits=[effmass.MassFit(carrier='hueco', band=0,
        direction='transversal 1', mass=-.321, r2=.9999, npts=9, window=.06)])
    original = copy.deepcopy(run)
    text = effmass.report(run)
    assert 'hole' in text and 'transverse 1' in text and '-0.321' in text
    assert 'hueco' not in text and 'transversal 1' not in text
    assert run == original
    table = Path(effmass.export(run, tmp_path)[0]).read_text(encoding='utf-8')
    row = next(line for line in table.splitlines() if not line.startswith('#'))
    assert row.split()[0] == 'hueco' and row.endswith('transversal 1')
    assert float(row.split()[2]) == -.321


def test_project_pdf_source_is_english_and_keeps_quality_codes(tmp_path):
    root, data = project.init(tmp_path / 'project', name='My sample')
    before = copy.deepcopy(data)
    gate = quality.evaluate(root, data)
    assert gate['verdict'] == 'bloqueado'
    text = '\n'.join(report._lines(root, data))
    assert 'reproducible project report' in text and 'BLOCKED' in text
    assert 'My sample' in text and 'Fuentes' not in text
    assert data == before


@pytest.mark.parametrize('kind', ['bands', 'combined'])
@pytest.mark.parametrize('ref,expected,shift', [
    ('fermi', r'$E - E_\mathrm{F}$ (eV)', -.8),
    ('vbm', r'$E - E_\mathrm{VBM}$ (eV)', -.8),
    ('none', r'$E$ (eV)', 0.),
])
def test_band_figures_label_the_actual_energy_reference(kind, ref, expected, shift, monkeypatch):
    from tests.test_benchmark_regressions import structure
    from qekit.modules import bands, combined, dos
    from qekit.core import style
    bs = structure(empty=True)
    bs.kdist = np.array([0., 1.])
    captured = []
    def save(figure, *args, **kwargs):
        axis = figure.axes[0]
        captured.append(axis.get_ylabel())
        assert axis.lines[0].get_ydata()[0] == pytest.approx(-2. - shift)
        return []
    monkeypatch.setattr(style, 'save', save)
    if kind == 'bands':
        bands.plot(bs, ref=ref, mark_extrema=False)
    else:
        dd = dos.DOSData(energies=np.linspace(-3., 3., 7), total=np.ones((1, 7)), fermi=-.8)
        combined.plot(bs, dd, ref=ref, mark_extrema=False)
    assert captured == [expected]
