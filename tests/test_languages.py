"""Unified language selection must preserve scripting and scientific settings."""
import json
import sys
from pathlib import Path

import pytest

from qekit import cli, config
from qekit.core import i18n
from qekit.modules import update


@pytest.fixture
def preferences(tmp_path, monkeypatch):
    monkeypatch.setattr(config, 'CONFIG_DIR', tmp_path)
    monkeypatch.setattr(config, 'CONFIG_FILE', tmp_path / 'config.ini')
    monkeypatch.setattr(config, '_migrar_si_hace_falta', lambda: False)
    monkeypatch.delenv('OLLA_DFT_LANG', raising=False)
    i18n.set_language(None)
    return config.CONFIG_FILE


def terminal(monkeypatch, answers):
    monkeypatch.setattr(sys.stdin, 'isatty', lambda: True)
    monkeypatch.setattr(sys.stdout, 'isatty', lambda: True)
    sequence = iter(answers)
    monkeypatch.setattr('builtins.input', lambda prompt: next(sequence))


@pytest.mark.parametrize('choice,language', [('1', 'en'), ('2', 'es'), ('', 'en')])
def test_startup_selection_and_persistence(preferences, monkeypatch, capsys, choice, language):
    config.set_value('ecutwfc', '85')
    config.set_value('pseudo_dir', '/custom/pseudos')
    terminal(monkeypatch, [choice, '0'])
    assert cli.main([]) == 0
    assert config.load()['language'] == language
    assert config.load()['ecutwfc'] == '85'
    assert config.load()['pseudo_dir'] == '/custom/pseudos'
    assert 'Language / Idioma' in capsys.readouterr().out
    i18n.set_language(None)
    assert i18n.get_language() == language


def test_invalid_selection_retries_and_saved_default_is_used(preferences, monkeypatch):
    config.set_value('language', 'es')
    terminal(monkeypatch, ['fr', '', '0'])
    assert cli.main([]) == 0
    assert config.load()['language'] == 'es'


def test_switch_language_inside_menu_overrides_environment(preferences, monkeypatch, capsys):
    monkeypatch.setenv('OLLA_DFT_LANG', 'en')
    terminal(monkeypatch, ['l', '2', 'p', '0'])
    seen = []
    monkeypatch.setitem(cli._DISPATCH, 'start', lambda args: seen.append(args.language) or 0)
    assert cli.main(['--language', 'en']) == 0
    assert seen == ['es']
    assert config.load()['language'] == 'es'
    assert cli._menu_labels('es')['goodbye'] in capsys.readouterr().out


def test_menu_catalog_uses_selected_language(preferences, monkeypatch, capsys):
    monkeypatch.setenv('OLLA_DFT_LANG', 'en')
    terminal(monkeypatch, ['2', 'c', 'eos', '0'])
    assert cli.main([]) == 0
    out = capsys.readouterr().out
    assert i18n.ui('options', 'es') + ':' in out
    assert i18n.ui('options', 'en') + ':' not in out


@pytest.mark.parametrize('argv', [['--language', 'es'], ['--language=en']])
def test_explicit_flag_skips_selector_without_saving(preferences, monkeypatch, capsys, argv):
    terminal(monkeypatch, ['0'])
    assert cli.main(argv) == 0
    assert not preferences.exists()
    assert 'Olla-DFT — Language / Idioma' not in capsys.readouterr().out


@pytest.mark.parametrize('stdin_tty,stdout_tty', [(False, True), (True, False), (False, False)])
def test_nonterminal_prints_help_without_reading_input(preferences, monkeypatch, capsys, stdin_tty, stdout_tty):
    monkeypatch.setattr(sys.stdin, 'isatty', lambda: stdin_tty)
    monkeypatch.setattr(sys.stdout, 'isatty', lambda: stdout_tty)
    monkeypatch.setattr('builtins.input', lambda _: pytest.fail('unexpected prompt'))
    assert cli.main([]) == 0
    assert '--help' in capsys.readouterr().out
    assert not preferences.exists()


@pytest.mark.parametrize('argv', [['--help'], ['--version'], ['eos', '--help', '--language', 'es']])
def test_direct_help_never_prompts_or_saves(preferences, monkeypatch, argv):
    terminal(monkeypatch, [])
    with pytest.raises(SystemExit) as exc:
        cli.main(argv)
    assert exc.value.code == 0
    assert not preferences.exists()


def test_language_priority_for_direct_commands(preferences, monkeypatch):
    config.set_value('language', 'es')
    seen = []
    monkeypatch.setitem(cli._DISPATCH, 'start', lambda args: seen.append(args.language) or 0)
    monkeypatch.setattr('builtins.input', lambda _: pytest.fail('unexpected prompt'))
    cli.main(['start'])
    monkeypatch.setenv('OLLA_DFT_LANG', 'en')
    cli.main(['start'])
    cli.main(['start', '--language', 'es'])
    cli.main(['start'])
    assert seen == ['es', 'en', 'es', 'en']


@pytest.mark.parametrize('failure', [OSError('read only'), config.configparser.Error('bad config')])
def test_unsavable_preference_still_opens_menu(preferences, monkeypatch, capsys, failure):
    def fail(*args):
        raise failure
    monkeypatch.setattr(config, 'set_value', fail)
    terminal(monkeypatch, ['2', '0'])
    assert cli.main([]) == 0
    out = capsys.readouterr()
    assert 'No se pudo guardar' in out.err
    assert cli._menu_labels('es')['goodbye'] in out.out


@pytest.mark.parametrize('exception', [EOFError, KeyboardInterrupt])
def test_cancel_selector_has_no_config_side_effect(preferences, monkeypatch, exception):
    terminal(monkeypatch, [])
    def cancel(prompt):
        raise exception
    monkeypatch.setattr('builtins.input', cancel)
    with pytest.raises(SystemExit) as exc:
        cli.main([])
    assert exc.value.code == 0
    assert not preferences.exists()


def test_failed_atomic_save_preserves_existing_config(preferences, monkeypatch):
    config.set_value('pseudo_dir', '/existing/pseudos')
    original = preferences.read_bytes()
    def fail(*args):
        raise OSError('no space')
    monkeypatch.setattr(config.os, 'replace', fail)
    with pytest.raises(OSError):
        config.set_value('language', 'es')
    assert preferences.read_bytes() == original
    assert list(preferences.parent.iterdir()) == [preferences]


@pytest.mark.parametrize('language', ['en', 'es'])
def test_one_update_source_for_both_languages(preferences, monkeypatch, language):
    i18n.set_language(language)
    monkeypatch.setattr(i18n, 'DEFAULT_LANGUAGE', language)
    monkeypatch.setattr(update, 'install_source', lambda: ('unknown', ''))
    plan = update.make_plan('1.4.0')
    assert 'olla-dft.git@v1.4.0' in plan.commands[0][-1] or 'olla-dft@v1.4.0' in plan.commands[0][-1]
    assert 'olla-dft-esp' not in json.dumps(plan.commands)


def test_language_catalogs_and_documentation_are_packaged():
    root = Path(__file__).resolve().parents[1]
    for lang in ('en', 'es'):
        for name in ('menu', 'studio', 'onboarding', 'dashboard', 'cli'):
            assert json.loads((root / f'qekit/data/i18n/{name}_{lang}.json').read_text())
    for name in ('README.es.md', 'docs/COMANDOS.md', 'docs/TEORIA.md', 'docs/COMMANDS.md', 'docs/THEORY.md'):
        assert (root / name).is_file()


def test_legacy_settings_survive_language_selection(preferences, monkeypatch):
    # Same INI schema as 1.3.x, including keys unknown to this release.
    preferences.write_text('[qekit]\npseudo_dir = /datos/cálculos/pseudos\n'
                           'ecutwfc = 90\nnproc = 8\nfuture_option = keep me\n'
                           '[external]\nsetting = untouched\n', encoding='utf-8')
    terminal(monkeypatch, ['2', '0'])
    assert cli.main([]) == 0
    saved = config.configparser.ConfigParser()
    saved.read(preferences, encoding='utf-8')
    assert saved['qekit']['pseudo_dir'] == '/datos/cálculos/pseudos'
    assert saved['qekit']['ecutwfc'] == '90'
    assert saved['qekit']['nproc'] == '8'
    assert saved['qekit']['future_option'] == 'keep me'
    assert saved['external']['setting'] == 'untouched'
    assert saved['qekit']['language'] == 'es'


def test_menu_switch_rebuilds_catalog_help(preferences, monkeypatch, capsys):
    terminal(monkeypatch, ['c', 'eos', 'l', '2', 'c', 'eos', '0'])
    assert cli.main(['--language', 'en']) == 0
    output = capsys.readouterr().out
    assert i18n.ui('options', 'en') + ':' in output
    assert i18n.ui('options', 'es') + ':' in output


@pytest.mark.parametrize('value,expected', [('ES', 'es'), ('es_ES.UTF-8', 'en'), ('fr', 'en')])
def test_environment_language_contract(preferences, monkeypatch, value, expected):
    monkeypatch.setenv('OLLA_DFT_LANG', value)
    assert i18n.set_language(None) == expected


def test_invalid_flag_does_not_prompt_or_save(preferences, monkeypatch):
    terminal(monkeypatch, [])
    assert cli.main(['--language', 'fr']) == 2
    assert not preferences.exists()
