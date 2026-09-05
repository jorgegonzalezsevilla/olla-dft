"""German is a real interface locale, with unchanged scientific payloads."""
import json
from pathlib import Path
import re

import pytest

from qekit import cli, config
from qekit.core import i18n
from qekit.modules import dashboard, docs, onboarding, recipes, studio, theory, wizard

ROOT = Path(__file__).resolve().parents[1]


def shape(value):
    if isinstance(value, dict):
        return {key: shape(item) for key, item in value.items()}
    if isinstance(value, list):
        return [shape(item) for item in value]
    return type(value).__name__


@pytest.mark.parametrize('name', ['cli', 'menu', 'dashboard', 'docs', 'onboarding', 'studio', 'recipes', 'wizard'])
def test_german_catalog_has_all_keys(name):
    base = ROOT / 'qekit/data/i18n'
    en = json.loads((base / f'{name}_en.json').read_text(encoding='utf-8'))
    de = json.loads((base / f'{name}_de.json').read_text(encoding='utf-8'))
    # Keywords may include additional German search terms.
    for catalog in (en, de):
        catalog.pop('keywords', None)
    assert shape(en) == shape(de)
    def placeholders(value):
        if isinstance(value, dict):
            return {key: placeholders(item) for key, item in value.items()}
        if isinstance(value, list):
            return [placeholders(item) for item in value]
        return sorted(re.findall(r"\{[A-Za-z_][^{}]*\}", value)) if isinstance(value, str) else []
    assert placeholders(en) == placeholders(de)


def test_german_language_selection_and_nested_command(tmp_path, monkeypatch, capsys):
    import sys
    monkeypatch.setattr(config, 'CONFIG_DIR', tmp_path)
    monkeypatch.setattr(config, 'CONFIG_FILE', tmp_path / 'config.ini')
    monkeypatch.setattr(config, '_migrar_si_hace_falta', lambda: False)
    monkeypatch.setenv('OLLA_DFT_LANG', 'en')
    monkeypatch.setattr(sys.stdin, 'isatty', lambda: True)
    monkeypatch.setattr(sys.stdout, 'isatty', lambda: True)
    answers = iter(['3', 'p', '0'])
    monkeypatch.setattr('builtins.input', lambda _: next(answers))
    seen = []
    monkeypatch.setitem(cli._DISPATCH, 'start', lambda args: seen.append(args.language) or 0)
    assert cli.main([]) == 0
    assert seen == ['de']
    assert config.load()['language'] == 'de'
    assert cli._menu_labels('de')['goodbye'] in capsys.readouterr().out


def test_german_help_and_sources_have_same_commands(capsys):
    with pytest.raises(SystemExit) as exc:
        cli.main(['eos', '--help', '--language', 'de'])
    assert exc.value.code == 0
    text = capsys.readouterr().out
    assert '--equation' in text and 'Optionen' in text
    assert {x.nombre for x in docs.extraer('en')} == {x.nombre for x in docs.extraer('de')}


def test_german_guides_preserve_commands_and_ids():
    for source, translated in zip(recipes.recetas('en'), recipes.recetas('de')):
        assert source.clave == translated.clave
        assert [p.comando for p in source.pasos] == [p.comando for p in translated.pasos]
        def paths(steps):
            return [[value for value in step.escribe if not any(c.isspace() for c in value) and re.search(r"[/]|\.[A-Za-z]", value)]
                    for step in steps]
        assert paths(source.pasos) == paths(translated.pasos)
    assert [m.clave for m in wizard.metas('en')] == [m.clave for m in wizard.metas('de')]
    assert wizard.buscar('Licht absorbieren', language='de')
    assert recipes.buscar('Bandlücke', language='de')
    assert onboarding._labels('de')['title'] != onboarding._labels('es')['title']


def test_german_explorer_keeps_original_scientific_values(tmp_path):
    rows = [{'id': 'sample', 'formula': 'Si', 'metrics': {
        'energy_per_atom': {'value': -123.456789012345, 'unit': 'eV/atom'}}}]
    payloads = []
    for language in ('en', 'es', 'de'):
        path = studio.generate(rows, tmp_path / f'{language}.html', language=language)
        text = path.read_text(encoding='utf-8')
        payload = json.loads(re.search(r'<script id="studio-data" type="application/json">(.*?)</script>', text, re.S).group(1))
        assert payload['language'] == language
        assert set(payload['labels']) == {'en', 'es', 'de'}
        assert '<option value="de">Deutsch</option>' in text
        payloads.append(payload['rows'])
    assert payloads[0] == payloads[1] == payloads[2]
    assert payloads[0][0]['metrics']['energy_per_atom']['value'] == -123.456789012345


def test_german_html_reference_has_localized_help_and_english_science(tmp_path):
    path = Path(docs.generar(str(tmp_path / 'reference.html'), language='de'))
    content = path.read_text(encoding='utf-8')
    assert '<html lang="de">' in content
    assert 'lang="en"' in content
    assert 'Birch' in content
    assert 'No hay comandos' not in content
    assert 'no hay fundamento' not in content
    assert len(theory.secciones('de')) == len(theory.secciones('en'))


def test_single_dashboard_does_not_link_missing_language_files(tmp_path):
    from qekit.modules import project
    root, data = project.init(tmp_path / 'project', name='Language test')
    path = dashboard.generate(root, data, language='de')
    assert 'class="language-switch"' not in path.read_text(encoding='utf-8')
    pages = dashboard.generate_all(root, data)
    assert len(pages) == 3
    for page in pages:
        content = page.read_text(encoding='utf-8')
        links = re.findall(r'class="language-switch"[^>]*href="([^"]+)"', content)
        assert len(links) == 2
        assert all((page.parent / name).is_file() for name in links)


@pytest.mark.parametrize("language,answer,expected", [("de", "ja", True), ("de", "nein", False), ("en", "yes", True), ("es", "sí", True)])
def test_yes_no_uses_interface_language(language, answer, expected, monkeypatch):
    monkeypatch.setattr(i18n, "get_language", lambda: language)
    prompts = []
    monkeypatch.setattr("builtins.input", lambda prompt: prompts.append(prompt) or answer)
    assert cli._ask_yes("Continue?") is expected
    assert {"de": "j/N", "en": "y/N", "es": "s/N"}[language] in prompts[0]
