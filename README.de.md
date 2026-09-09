<h1 align="center">Olla-DFT</h1>

<p align="center"><b>Von Kristallstrukturen zu Ergebnissen mit Quantum ESPRESSO.</b><br>
Berechnungen vorbereiten, Materialeigenschaften auswerten und Abbildungen teilen.</p>

<p align="center"><a href="README.md">English</a> · <a href="README.es.md">Español</a> · <b>Deutsch</b></p>

<p align="center"><img src="examples/demo_Si/Si_bandas_dos.png" width="820" alt="Elektronische Bandstruktur und Zustandsdichte von Silizium"><br><sub>Berechnetes Beispiel mit Quantum ESPRESSO. Die LDA-Bandlücke entspricht nicht der experimentellen Bandlücke.</sub></p>

## Möglichkeiten

| Bereich | Funktionen |
|---|---|
| Vorbereitung | Strukturen, Symmetrie, Pseudopotentiale, k-Punkt-Netze und Bandpfade. |
| Elektronen und Spektren | Bänder, DOS/PDOS, Bandlücken, Magnetismus und optische Eigenschaften. |
| Materialeigenschaften | Phononen, Thermodynamik, Zustandsgleichungen, Elastizität, Oberflächen und Defekte. |
| Projekte | Geführter Einstieg, Berechnungsserien, Qualitätsprüfungen und nachvollziehbare Ergebnisse. |
| Darstellung | Anpassbare Abbildungen, interaktiver Offline-Explorer und Datenexport. |

[Alle Befehle](https://github.com/jorgegonzalezsevilla/olla-dft/blob/main/docs/BEFEHLE.md) · [Wissenschaftliche Grundlagen (Englisch)](https://github.com/jorgegonzalezsevilla/olla-dft/blob/main/docs/THEORY.md) · [Validierung (Englisch)](https://github.com/jorgegonzalezsevilla/olla-dft/blob/main/docs/VALIDATION.md)

## Installation und Sprache

Erfordert **Python 3.9+**. Für neue Berechnungen sind Quantum ESPRESSO und
Pseudopotentiale zusätzlich erforderlich. Vorhandene Ergebnisse lassen sich
auch ohne Quantum ESPRESSO auswerten.

```bash
git clone https://github.com/jorgegonzalezsevilla/olla-dft.git
cd olla-dft
python3 -m venv .venv
source .venv/bin/activate
pip install .
olla-dft
```

Wähle **3 — Deutsch**. Die Auswahl wird gespeichert; mit **l** lässt sich die
Sprache im Menü ändern. `olla-dft --language de` öffnet das Menü direkt auf
Deutsch. Für Skripte funktionieren `--language de`, `OLLA_DFT_LANG=de` oder
`olla-dft config set language de`, ohne zusätzliche Rückfragen.

**Wissenschaftliche Berichte werden auf Englisch erstellt**, unabhängig von
der Oberflächensprache. Die ausführliche Theorie wird für die deutsche
Oberfläche ebenfalls auf Englisch bereitgestellt. Befehle, Dateiformate und
wissenschaftliche Daten behalten ihre ursprünglichen Bezeichner.

## Ergebnisse untersuchen und exportieren

```bash
olla-dft results ingest ./calculation --project ./my-project
olla-dft results explore --project ./my-project --language de -o results.html
```

Öffne `results.html`, filtere Berechnungen, wähle Größen und Einheiten und passe
Titel, Farbe und Achsen an. Exportiere **SVG, PNG, CSV, JSON oder interaktives
HTML**. Der Explorer arbeitet offline mit einer festen Kopie der Daten;
erzeuge ihn erneut, wenn weitere Ergebnisse hinzugekommen sind.

[Beispiele](https://github.com/jorgegonzalezsevilla/olla-dft/tree/main/examples/) · [Galerie](https://jorgegonzalezsevilla.github.io/olla-dft-bench/publication-1.2.0/index-en.html) · [Exporthinweise (Englisch)](https://github.com/jorgegonzalezsevilla/olla-dft/blob/main/docs/RESULTS-EXPLORER.md)

## Nach einer Unterbrechung weiterrechnen

`olla-dft resilient` schreibt und prüft Checkpoints, um unterstützte
`pw.x`-Rechnungen fortzusetzen, solange die Festplatte erhalten bleibt. Dafür
muss zuerst die persistente Umgebung eingerichtet werden:
[Wiederherstellungsanleitung (Spanisch)](https://github.com/jorgegonzalezsevilla/olla-dft/blob/main/docs/resilience/RECUPERACION.md).

Geprüft wurden lokale SCF-, `relax`- und `vc-relax`-Paare mit simulierten
Prozessabbrüchen. Die Wiederherstellung nach einem physischen Stromausfall oder
einem Festplattenverlust **ist nicht nachgewiesen**.
[Ergebnisse und Toleranzen](https://jorgegonzalezsevilla.github.io/olla-dft-bench/publication-1.2.0/index-en.html) ·
[Wiederherstellungsvertrag (Englisch)](https://github.com/jorgegonzalezsevilla/olla-dft/blob/main/docs/resilience/CONTRACT.md).

## Dokumentation und Lizenz

Ein Projekt von **Jorge Enrique González Sevilla**, unabhängig von Quantum
ESPRESSO. Freie Software unter **AGPL-3.0-or-later** mit den dokumentierten
Ausnahmen für einzelne Dateien. Keine automatische Telemetrie.

[Befehle](https://github.com/jorgegonzalezsevilla/olla-dft/blob/main/docs/BEFEHLE.md) · [Wissenschaftliche Referenz (Englisch)](https://github.com/jorgegonzalezsevilla/olla-dft/blob/main/docs/THEORY.md) · [Validierung (Englisch)](https://github.com/jorgegonzalezsevilla/olla-dft/blob/main/docs/VALIDATION.md) ·
[Sprachen und Migration](https://github.com/jorgegonzalezsevilla/olla-dft/blob/main/docs/LANGUAGES.md) · [Änderungen](https://github.com/jorgegonzalezsevilla/olla-dft/blob/main/CHANGELOG.md) ·
[Fehler melden](https://github.com/jorgegonzalezsevilla/olla-dft/issues) ·
[Lizenz](https://github.com/jorgegonzalezsevilla/olla-dft/blob/main/LICENSE) · [Lizenzumfang](https://github.com/jorgegonzalezsevilla/olla-dft/blob/main/LICENSING.md)

Zitiere die verwendete Version über [CITATION.cff](https://github.com/jorgegonzalezsevilla/olla-dft/blob/main/CITATION.cff) und den
[einen Softwareeintrag auf Zenodo](https://doi.org/10.5281/zenodo.22263121).
Zitiere außerdem Quantum ESPRESSO und die verwendeten Pseudopotentiale.
