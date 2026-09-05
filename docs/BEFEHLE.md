# Olla-DFT-Befehlsreferenz

Die 80 Unterbefehle von `olla-dft`, nach Themen gruppiert. Aus dem Programmcode erzeugt mit `python tools/build_docs.py`. Hilfe: `olla-dft COMMAND --help --language de`; HTML-Referenz: `olla-dft docs --language de`. Wissenschaftliche Berichte und die ausführliche Theorie sind auf Englisch.

## Inhalt

- **Erste Schritte**: [`start`](#start), [`wizard`](#wizard), [`recetas`](#recetas), [`teoria`](#teoria), [`docs`](#docs), [`sistema`](#sistema), [`selftest`](#selftest), [`update`](#update)
- **Strukturen und Eingaben**: [`gen`](#gen), [`info`](#info), [`kpath`](#kpath), [`prim`](#prim), [`conv`](#conv), [`supercell`](#supercell), [`convert`](#convert)
- **Elektronische Struktur**: [`bands`](#bands), [`dos`](#dos), [`plot`](#plot), [`gap`](#gap), [`fermi`](#fermi), [`effmass`](#effmass), [`wannier`](#wannier), [`unfold`](#unfold), [`topology`](#topology), [`hubbard`](#hubbard)
- **Spektren und Antwort**: [`optics`](#optics), [`tddft`](#tddft), [`xanes`](#xanes), [`xps`](#xps), [`corehole`](#corehole), [`charge`](#charge), [`charges`](#charges), [`wf`](#wf), [`berry`](#berry)
- **Phononen, Transport und Temperatur**: [`phonons`](#phonons), [`elph`](#elph), [`transport`](#transport), [`ballistic`](#ballistic), [`kappa`](#kappa), [`qha`](#qha), [`thermochem`](#thermochem), [`md`](#md), [`derived`](#derived)
- **Mechanik und Stabilität**: [`converge`](#converge), [`eos`](#eos), [`elastic`](#elastic), [`strain`](#strain), [`layers`](#layers), [`xrd`](#xrd), [`exfoliate`](#exfoliate), [`gamma`](#gamma)
- **Oberflächen, Defekte und Chemie**: [`surface`](#surface), [`defect`](#defect), [`interface`](#interface), [`adsorb`](#adsorb), [`eform`](#eform), [`align`](#align), [`esm`](#esm), [`echem`](#echem), [`neb`](#neb), [`amorphous`](#amorphous)
- **Automatisierung und Qualität**: [`doctor`](#doctor), [`audit`](#audit), [`crosscheck`](#crosscheck), [`cost`](#cost), [`db`](#db), [`hull`](#hull), [`mlip`](#mlip), [`suggest`](#suggest), [`datasheet`](#datasheet), [`report`](#report), [`compare`](#compare), [`tune`](#tune), [`results`](#results), [`campaign`](#campaign), [`pseudos`](#pseudos)
- **Projekt**: [`project`](#project), [`resilient`](#resilient)
- **Erscheinungsbild und Konfiguration**: [`templates`](#templates), [`config`](#config)

## Erste Schritte

### `start`

geführter Start, um ein Projekt ohne CLI-Kenntnisse anzulegen

**Verwendung:** `olla-dft start [-h] [--project PROJECT] [--structure STRUCTURE] [--goal GOAL] [--name NAME] [--non-interactive] [--no-validate] [--language {es,en,de}]`

**Optionen:**

| Option | Beschreibung |
|---|---|
| `--project` | Projektordner (Standard: `.`) |
| `--structure` | CIF, POSCAR oder pw.x-Eingabe |
| `--goal` | relax, gap, dos, phonons, optics oder scf |
| `--name` | sichtbarer Name des Projekts |
| `--non-interactive` | nicht nachfragen; erfordert --structure in einem neuen Projekt |
| `--no-validate` | die anfängliche Validierung nicht ausführen |
| `--language {es,en,de}` | Oberflächensprache (en, es, de; gemäß --language) |

### `wizard`

Assistent: sag mir, WAS du wissen willst, und ich sage dir, was zu rechnen ist, in Reihenfolge und mit den Befehlen

**Verwendung:** `olla-dft wizard [-h] [--goal GOAL] [--ask TEXTO] [--list] [--term TERM] [--no-glossary] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [file]`

**Argumente:**

- `file` — deine Struktur (optional)

**Optionen:**

| Option | Beschreibung |
|---|---|
| `--goal` | Schlüssel des Ziels; mit --list werden sie aufgelistet |
| `--ask TEXTO` | beschreibe dein Ziel auf Deutsch, Englisch oder Spanisch, z. B. 'Licht absorbieren' |
| `--list` | alles auflisten, was der Assistent kann |
| `--term` | was ein Begriff bedeutet |
| `--no-glossary` | die Fachbegriffe am Ende der Antwort nicht erklären |
| `--pseudo-dir` | Ordner mit den UPF-Pseudopotentialen (ohne Angabe der aus 'olla-dft config') |
| `--pseudo EL=UPF` | ein bestimmtes Pseudopotential erzwingen, z. B. Fe=Fe.rel-pbe.UPF. Kann wiederholt werden. Ohne diese Angabe wählt Olla-DFT mit 'olla-dft pseudos' |

### `recetas`

vollständige Sitzungen von Anfang bis Ende: welcher Befehl nach welchem kommt und welche Datei sie sich übergeben

**Verwendung:** `olla-dft recetas [-h] [--buscar TEXTO] [--script [ARCHIVO]] [receta]`

**Argumente:**

- `receta` — Schlüssel des Rezepts; ohne Angabe werden alle aufgelistet

**Optionen:**

| Option | Beschreibung |
|---|---|
| `--buscar TEXTO` | mit eigenen Worten suchen, ohne den Schlüssel zu kennen |
| `--script ARCHIVO` | das Rezept als kommentiertes Shell-Skript schreiben, bereit zum Bearbeiten |

### `teoria`

die physikalische Grundlage eines Befehls: was er beantwortet, die Formeln, die er implementiert, aus welchem Modul sie stammen und woher jede Zahl kommt

**Verwendung:** `olla-dft teoria [-h] [--all] [-o ARCHIVO.md] [comando]`

**Argumente:**

- `comando` — zu erklärender Befehl; ohne Angabe der Index

**Optionen:**

| Option | Beschreibung |
|---|---|
| `--all` | das vollständige Dokument (alle Bereiche) |
| `-o, --output ARCHIVO.md` | als Markdown speichern, statt es auszugeben |

### `docs`

navigierbare Referenz aller Unterbefehle, aus dem Code selbst erzeugt

**Verwendung:** `olla-dft docs [-h] [-o OUTPUT] [--open] [--language {es,en,de}] [--both] [--all-languages]`

**Optionen:**

| Option | Beschreibung |
|---|---|
| `-o, --output` | HTML-Ausgabedatei (Standard: `olla-dft-docs.html`) |
| `--open` | nach Abschluss im Browser öffnen |
| `--language {es,en,de}` | Oberflächensprache (en, es, de; gemäß --language) |
| `--both` | getrennte Referenzen auf Spanisch und Englisch erzeugen |
| `--all-languages` | Oberflächen auf Englisch, Spanisch und Deutsch erzeugen |

### `sistema`

was Olla-DFT auf dieser Maschine sieht: Kodierung, wo es die Konfiguration speichert, welche QE-Binärdateien es findet und wie Rechnungen hier gestartet werden

**Verwendung:** `olla-dft sistema [-h]`

### `selftest`

Olla-DFT gegen veröffentlichte Werte prüfen, nicht gegen sich selbst

**Verwendung:** `olla-dft selftest [-h] [--full] [--mlip] [--only ONLY] [--list] [--pseudo-dir PSEUDO_DIR] [--pw-cmd PW_CMD] [--nproc NPROC] [-j JOBS] [--keep CARPETA]`

**Optionen:**

| Option | Beschreibung |
|---|---|
| `--full` | die Tests einschließen, die pw.x wirklich ausführen (etwa zehn Minuten) |
| `--mlip` | den Test mit gelerntem Potential separat einschließen (erfordert MACE) |
| `--only` | nur diese Tests, kommagetrennt |
| `--list` | die Tests und ihre Referenzen auflisten, ohne etwas auszuführen |
| `--pseudo-dir` | Pseudopotentiale für die --full-Tests |
| `--pw-cmd` | pw.x-Programm für --run; aus seinem Pfad werden die übrigen QE-Binärdateien gefunden |
| `--nproc NPROC` | Anzahl der MPI-Prozesse für die mit --run gestarteten Rechnungen |
| `-j, --jobs JOBS` | gleichzeitige Tests (Standard: 1) |
| `--keep CARPETA` | die Rechnungen hier behalten, statt sie zu löschen |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria selftest`](THEORY.md)

### `update`

prüfen, ob eine neuere Version von Olla-DFT veröffentlicht ist, und sie gegebenenfalls nach einer Bestätigung installieren; läuft nie von allein

**Verwendung:** `olla-dft update [-h] [--check] [--yes] [--version TAG]`

**Optionen:**

| Option | Beschreibung |
|---|---|
| `--check` | nur prüfen und berichten, nichts installieren |
| `--yes` | nicht nachfragen; direkt installieren, wenn eine neuere Version existiert |
| `--version TAG` | eine bestimmte Version (z. B. v1.0.1) statt der neuesten installieren |

## Strukturen und Eingaben

### `gen`

pw.x-Eingaben und Post-Processing erzeugen

**Verwendung:** `olla-dft gen [-h] [-p {scf,relax,vc-relax,nscf,bands,dos,all,md}] [-o OUTDIR] [-k {coarse,fine,gamma,medium,very-fine}] [--kspacing KSPACING] [--kgrid N N N] [--band-points BAND_POINTS] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--insulator] [--primitive] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--prefix PREFIX] [--nspin {1,2}] [--mag MAG] [--vdw {grimme-d2,grimme-d3,DFT-D,ts-vdw,xdm,mbd}] [--soc] [--hubbard EL=U] [--hubbard-style {legacy,card}] [--charge Q] [--dipole [EJE]] [--nosym] [--functional {b3lyp,gaupbe,hse,pbe0}] [--exx-grid NxNxN] [--exx-fraction EXX_FRACTION] [--dt FS] [--nstep NSTEP] [--thermostat {none,rescaling,berendsen,andersen,initial,reduce-history}] [-T TEMPERATURE] file`

**Argumente:**

- `file` — Struktur (CIF, POSCAR, pw.x-Eingabe, ...)

**Optionen:**

| Option | Beschreibung |
|---|---|
| `-p, --preset {scf,relax,vc-relax,nscf,bands,dos,all,md}` | Rechnungstyp (Standard: scf) |
| `-o, --outdir` | Ausgabeordner (Standard: `.`) |
| `-k, --klevel {coarse,fine,gamma,medium,very-fine}` | Dichte des k-Gitters (gamma/coarse/medium/fine/very-fine) |
| `--kspacing KSPACING` | k-Abstand in Å^-1 (überschreibt --klevel) |
| `--kgrid N` | explizites k-Gitter für scf/relax (drei ganze Zahlen; überschreibt --kspacing und --klevel) |
| `--band-points BAND_POINTS` | Punkte pro Segment des k-Pfads |
| `--ecutwfc ECUTWFC` | Cutoff der Wellenfunktionen (Ry) |
| `--ecutrho ECUTRHO` | Cutoff der Dichte (Ry) |
| `--insulator` | occupations='fixed' (Isolatoren; Standard: smearing) |
| `--primitive` | vor dem Erzeugen auf die standardisierte primitive Zelle reduzieren |
| `--pseudo-dir` | Pseudopotential-Ordner (überschreibt config) |
| `--pseudo EL=UPF` | ein bestimmtes Pseudopotential erzwingen, z. B. Fe=Fe.rel-pbe.UPF. Kann wiederholt werden. Ohne diese Angabe wählt Olla-DFT mit 'olla-dft pseudos' |
| `--prefix` | Prefix der Rechnung (Standard: Formel) |
| `--nspin {1,2}` | 2 schaltet die Spinpolarisation ein (Standard: `1`) |
| `--mag` | Startmagnetisierung: eine Zahl (0.5) oder pro Element (Fe=0.7,O=0). Impliziert --nspin 2 |
| `--vdw {grimme-d2,grimme-d3,DFT-D,ts-vdw,xdm,mbd}` | Dispersionskorrektur (van der Waals) |
| `--soc` | Spin-Bahn-Kopplung: nichtkollineare Rechnung mit lspinorb (erfordert vollrelativistische Pseudopotentiale) |
| `--hubbard EL=U` | Hubbard-U in eV pro Element, z. B. Ni=4.1. Kann wiederholt werden. Um es zu berechnen, statt es zu raten:  olla-dft hubbard --cycle |
| `--hubbard-style {legacy,card}` | legacy = lda_plus_u (QE <= 7.0), card = HUBBARD-Karte (QE >= 7.1) (Standard: `legacy`) |
| `--charge Q` | Gesamtladung der Zelle (tot_charge): +1 entfernt ein Elektron, -1 fügt eines hinzu |
| `--dipole EJE` | Dipolkorrektur für polare Slabs; ohne Wert wird die c-Achse verwendet. Legt den Sägezahn ins Vakuum |
| `--nosym` | Symmetrie ausschalten (nosym und noinv) |
| `--functional {b3lyp,gaupbe,hse,pbe0}` | Hybridfunktional: hse, pbe0, b3lyp oder gaupbe. Kostet ein bis zwei Größenordnungen mehr als PBE, und der Bericht sagt das mit Zahlen |
| `--exx-grid NxNxN` | q-Gitter für den exakten Austausch (Standard 1x1x1). Muss das k-Gitter teilen |
| `--exx-fraction EXX_FRACTION` | Anteil des exakten Austauschs, falls du den des Funktionals ändern willst |
| `--dt FS` | MD-Zeitschritt in fs (Standard: 1.0) |
| `--nstep NSTEP` | MD-Schritte (Standard: 1000) |
| `--thermostat {none,rescaling,berendsen,andersen,initial,reduce-history}` | MD-Thermostat; none = NVE (Standard) |
| `-T, --temperature TEMPERATURE` | Zieltemperatur der MD in K (Standard: 300) |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria gen`](THEORY.md)

### `info`

Struktur- und Symmetrieinformationen

**Verwendung:** `olla-dft info [-h] file`

**Argumente:**

- `file` — Eingabestruktur (CIF, POSCAR, pw.x-Eingabe...)

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria info`](THEORY.md)

### `kpath`

Hochsymmetriepfad (seekpath)

**Verwendung:** `olla-dft kpath [-h] file`

**Argumente:**

- `file` — Eingabestruktur (CIF, POSCAR, pw.x-Eingabe...)

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria kpath`](THEORY.md)

### `prim`

standardisierte primitive Zelle

**Verwendung:** `olla-dft prim [-h] [-o OUTPUT] file`

**Argumente:**

- `file` — Eingabestruktur (CIF, POSCAR, pw.x-Eingabe...)

**Optionen:**

| Option | Beschreibung |
|---|---|
| `-o, --output` | Ausgabestrukturdatei (standardmäßig eine .cif mit dem Namen des Befehls) (Standard: `primitive.cif`) |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria prim`](THEORY.md)

### `conv`

standardisierte konventionelle Zelle

**Verwendung:** `olla-dft conv [-h] [-o OUTPUT] file`

**Argumente:**

- `file` — Eingabestruktur (CIF, POSCAR, pw.x-Eingabe...)

**Optionen:**

| Option | Beschreibung |
|---|---|
| `-o, --output` | Ausgabestrukturdatei (standardmäßig eine .cif mit dem Namen des Befehls) (Standard: `conventional.cif`) |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria conv`](THEORY.md)

### `supercell`

eine Superzelle bauen

**Verwendung:** `olla-dft supercell [-h] [-o OUTPUT] file nx ny nz`

**Argumente:**

- `file` — Eingabestruktur (CIF, POSCAR, pw.x-Eingabe...)
- `nx` — Wiederholungen der Zelle entlang a
- `ny` — Wiederholungen der Zelle entlang b
- `nz` — Wiederholungen der Zelle entlang c

**Optionen:**

| Option | Beschreibung |
|---|---|
| `-o, --output` | Ausgabestrukturdatei (standardmäßig eine .cif mit dem Namen des Befehls) (Standard: `supercell.cif`) |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria supercell`](THEORY.md)

### `convert`

Format konvertieren (CIF/POSCAR/XYZ)

**Verwendung:** `olla-dft convert [-h] [-o OUTPUT_FLAG] file [output]`

**Argumente:**

- `file` — Eingabestruktur (CIF, POSCAR, pw.x-Eingabe...)
- `output` — Zieldatei; das Format wird aus der Endung abgeleitet (.cif, .vasp, .xyz...)

**Optionen:**

| Option | Beschreibung |
|---|---|
| `-o, --output-flag` | Ausgabedatei (Alternative zur positionalen Angabe) |

## Elektronische Struktur

### `bands`

die Bandstruktur analysieren und plotten

**Verwendung:** `olla-dft bands [-h] [-o OUTDIR] [--prefix PREFIX] [--ref {auto,fermi,vbm,none}] [--emin EMIN] [--emax EMAX] [--no-plot] [--dpi DPI] [--format FORMAT] [-t TEMPLATE] [--size {paper,poster,presentation}] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--aspect ASPECT] [--mono] [--dashes {auto,always,never}] [--title TITLE] [--gap-label] [--panel PANEL] [--fat SELECTOR] [--fat-scale FAT_SCALE] [--projwfc ARCHIVO] [path]`

**Argumente:**

- `path` — Ordner der Rechnung (oder Pfad zur .xml)

**Optionen:**

| Option | Beschreibung |
|---|---|
| `-o, --outdir` | Ausgabeordner (Standard: `.`) |
| `--prefix` | Prefix der Rechnung (wird automatisch erkannt) |
| `--ref {auto,fermi,vbm,none}` | Energienullpunkt (Standard: auto) |
| `--emin EMIN` | untere Grenze der Energieachse (eV) (Standard: `-6.0`) |
| `--emax EMAX` | obere Grenze der Energieachse (eV) (Standard: `6.0`) |
| `--no-plot` | nur die Daten exportieren, ohne die Grafik zu erzeugen |
| `--dpi DPI` | Auflösung der Bitmap-Formate (Standard: `600`) |
| `--format` | kommagetrennte Formate: pdf,png,svg,eps,tif (Standard: `pdf,png`) |
| `-t, --template` | visuelle Vorlage: dark, journal, latex, latex-true, minimal, mono, mono-latex, poster, slides (oder der Pfad zu einer eigenen JSON-Datei) |
| `--size {paper,poster,presentation}` | Schriftskala: paper / presentation / poster |
| `--font {sans,serif,latex}` | Schriftfamilie (latex = Computer Modern) |
| `--usetex` | den Text mit echtem LaTeX rendern |
| `--palette` | Palette: grayscale, okabe-ito, okabe-ito-dark oder kommagetrennte Hexfarben |
| `--background` | Hintergrundfarbe, z. B. '#FFFFFF' oder 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | Spaltenbreiten des Verlags (Standard: `generic`) |
| `--width` | Breite: single / onehalf / double oder eine Zahl in mm |
| `--aspect ASPECT` | Verhältnis Höhe/Breite der Abbildung |
| `--mono` | monochrom: schwarze Tinte und Linienmuster (für Zeitschriften, die Farbe berechnen) |
| `--dashes {auto,always,never}` | Linienmuster als sekundäre Kodierung (Standard: `auto`) |
| `--title` | Titel innerhalb der Abbildung (standardmäßig keiner: in einem Artikel gehört der Text in die Bildunterschrift) |
| `--gap-label` | den Wert der Bandlücke in der Grafik annotieren |
| `--panel` | Panel-Beschriftung, z. B. '(a)' |
| `--fat SELECTOR` | Fatbands: Gewicht eines Orbitals auf jedem Band. Zum Beispiel Ni-d, Si-p, O, d oder atomo:3. Benötigt die projwfc.x-Ausgabe DERSELBEN Bandrechnung |
| `--fat-scale FAT_SCALE` | Größe der Fatband-Marker (Standard: `55.0`) |
| `--projwfc ARCHIVO` | projwfc.x-Ausgabe (standardmäßig projwfc.out im selben Ordner) |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria bands`](THEORY.md)

### `dos`

DOS und PDOS analysieren und plotten

**Verwendung:** `olla-dft dos [-h] [--mode {orbital,element,total}] [-o OUTDIR] [--prefix PREFIX] [--ref {auto,fermi,vbm,none}] [--emin EMIN] [--emax EMAX] [--no-plot] [--dpi DPI] [--format FORMAT] [-t TEMPLATE] [--size {paper,poster,presentation}] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--aspect ASPECT] [--mono] [--dashes {auto,always,never}] [--title TITLE] [--gap-label] [--panel PANEL] [--dband EL[-ORB]] [--dband-emax eV] [path]`

**Argumente:**

- `path` — Ordner der Rechnung (oder Pfad zur .xml)

**Optionen:**

| Option | Beschreibung |
|---|---|
| `--mode {orbital,element,total}` | wie die PDOS zerlegt wird (Standard: `orbital`) |
| `-o, --outdir` | Ausgabeordner (Standard: `.`) |
| `--prefix` | Prefix der Rechnung (wird automatisch erkannt) |
| `--ref {auto,fermi,vbm,none}` | Energienullpunkt (Standard: auto) |
| `--emin EMIN` | untere Grenze der Energieachse (eV) (Standard: `-6.0`) |
| `--emax EMAX` | obere Grenze der Energieachse (eV) (Standard: `6.0`) |
| `--no-plot` | nur die Daten exportieren, ohne die Grafik zu erzeugen |
| `--dpi DPI` | Auflösung der Bitmap-Formate (Standard: `600`) |
| `--format` | kommagetrennte Formate: pdf,png,svg,eps,tif (Standard: `pdf,png`) |
| `-t, --template` | visuelle Vorlage: dark, journal, latex, latex-true, minimal, mono, mono-latex, poster, slides (oder der Pfad zu einer eigenen JSON-Datei) |
| `--size {paper,poster,presentation}` | Schriftskala: paper / presentation / poster |
| `--font {sans,serif,latex}` | Schriftfamilie (latex = Computer Modern) |
| `--usetex` | den Text mit echtem LaTeX rendern |
| `--palette` | Palette: grayscale, okabe-ito, okabe-ito-dark oder kommagetrennte Hexfarben |
| `--background` | Hintergrundfarbe, z. B. '#FFFFFF' oder 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | Spaltenbreiten des Verlags (Standard: `generic`) |
| `--width` | Breite: single / onehalf / double oder eine Zahl in mm |
| `--aspect ASPECT` | Verhältnis Höhe/Breite der Abbildung |
| `--mono` | monochrom: schwarze Tinte und Linienmuster (für Zeitschriften, die Farbe berechnen) |
| `--dashes {auto,always,never}` | Linienmuster als sekundäre Kodierung (Standard: `auto`) |
| `--title` | Titel innerhalb der Abbildung (standardmäßig keiner: in einem Artikel gehört der Text in die Bildunterschrift) |
| `--gap-label` | den Wert der Bandlücke in der Grafik annotieren |
| `--panel` | Panel-Beschriftung, z. B. '(a)' |
| `--dband EL[-ORB]` | Zentrum, Breite und Füllung eines projizierten Bandes, z. B. Pt (verwendet d) oder Ni-p. Das ist der Deskriptor, der mit der Adsorptionsenergie korreliert |
| `--dband-emax eV` | obere Grenze des Integrals, relativ zum Fermi-Niveau |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria dos`](THEORY.md)

### `plot`

kombinierter Plot Bänder + DOS

**Verwendung:** `olla-dft plot [-h] [--mode {orbital,element,total}] [-o OUTDIR] [--prefix PREFIX] [--ref {auto,fermi,vbm,none}] [--emin EMIN] [--emax EMAX] [--no-plot] [--dpi DPI] [--format FORMAT] [-t TEMPLATE] [--size {paper,poster,presentation}] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--aspect ASPECT] [--mono] [--dashes {auto,always,never}] [--title TITLE] [--gap-label] [--panel PANEL] [path]`

**Argumente:**

- `path` — Ordner der Rechnung (oder Pfad zur .xml)

**Optionen:**

| Option | Beschreibung |
|---|---|
| `--mode {orbital,element,total}` | wie die PDOS zerlegt wird (Standard: `orbital`) |
| `-o, --outdir` | Ausgabeordner (Standard: `.`) |
| `--prefix` | Prefix der Rechnung (wird automatisch erkannt) |
| `--ref {auto,fermi,vbm,none}` | Energienullpunkt (Standard: auto) |
| `--emin EMIN` | untere Grenze der Energieachse (eV) (Standard: `-6.0`) |
| `--emax EMAX` | obere Grenze der Energieachse (eV) (Standard: `6.0`) |
| `--no-plot` | nur die Daten exportieren, ohne die Grafik zu erzeugen |
| `--dpi DPI` | Auflösung der Bitmap-Formate (Standard: `600`) |
| `--format` | kommagetrennte Formate: pdf,png,svg,eps,tif (Standard: `pdf,png`) |
| `-t, --template` | visuelle Vorlage: dark, journal, latex, latex-true, minimal, mono, mono-latex, poster, slides (oder der Pfad zu einer eigenen JSON-Datei) |
| `--size {paper,poster,presentation}` | Schriftskala: paper / presentation / poster |
| `--font {sans,serif,latex}` | Schriftfamilie (latex = Computer Modern) |
| `--usetex` | den Text mit echtem LaTeX rendern |
| `--palette` | Palette: grayscale, okabe-ito, okabe-ito-dark oder kommagetrennte Hexfarben |
| `--background` | Hintergrundfarbe, z. B. '#FFFFFF' oder 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | Spaltenbreiten des Verlags (Standard: `generic`) |
| `--width` | Breite: single / onehalf / double oder eine Zahl in mm |
| `--aspect ASPECT` | Verhältnis Höhe/Breite der Abbildung |
| `--mono` | monochrom: schwarze Tinte und Linienmuster (für Zeitschriften, die Farbe berechnen) |
| `--dashes {auto,always,never}` | Linienmuster als sekundäre Kodierung (Standard: `auto`) |
| `--title` | Titel innerhalb der Abbildung (standardmäßig keiner: in einem Artikel gehört der Text in die Bildunterschrift) |
| `--gap-label` | den Wert der Bandlücke in der Grafik annotieren |
| `--panel` | Panel-Beschriftung, z. B. '(a)' |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria plot`](THEORY.md)

### `gap`

nur der Bandlücken-Bericht (schnell)

**Verwendung:** `olla-dft gap [-h] [--prefix PREFIX] [path]`

**Argumente:**

- `path` — Ordner der Rechnung (oder Pfad zur .xml)

**Optionen:**

| Option | Beschreibung |
|---|---|
| `--prefix` | Prefix der Rechnung (wird automatisch erkannt) |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria gap`](THEORY.md)

### `fermi`

die Fermi-Fläche als BXSF exportieren

**Verwendung:** `olla-dft fermi [-h] [-o OUTDIR]`

**Optionen:**

| Option | Beschreibung |
|---|---|
| `-o, --outdir` | Ausgabeordner (Standard: `transporte`) |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria fermi`](THEORY.md)

### `effmass`

effektive Masse durch parabolischen Fit der Bänder

**Verwendung:** `olla-dft effmass [-h] [-o OUTDIR] [--bands-dir BANDS_DIR] [--collect] [--run] [--half-width HALF_WIDTH] [--points POINTS] [--window WINDOW] [--min-points MIN_POINTS] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--pw-cmd PW_CMD] [--nproc NPROC] [--timeout TIMEOUT] file`

**Argumente:**

- `file` — Eingabestruktur (CIF, POSCAR, pw.x-Eingabe...)

**Optionen:**

| Option | Beschreibung |
|---|---|
| `-o, --outdir` | Ausgabeordner (Standard: `masa_efectiva`) |
| `--bands-dir` | Ordner mit einer bereits durchgeführten Bandrechnung (daraus kommen VBM und CBM) |
| `--collect` | die bereits gelaufene feine Rechnung einlesen |
| `--run` | die feine Rechnung ausführen, sobald sie vorbereitet ist |
| `--half-width HALF_WIDTH` | Halbbreite jeder Linie in Å⁻¹ (Standard: `0.06`) |
| `--points POINTS` | k-Punkte pro Linie (ungerade) (Standard: `21`) |
| `--window WINDOW` | Halbbreite des Schnellfits entlang des Pfads, in Å⁻¹ zu jeder Seite des Extremums (standardmäßig die Hälfte der parabolischen Grenze: ±0.06) |
| `--min-points MIN_POINTS` | Mindestpunkte für den Schnellfit (Standard: `7`) |
| `--pseudo-dir` | Ordner mit den UPF-Pseudopotentialen (ohne Angabe der aus 'olla-dft config') |
| `--pseudo EL=UPF` | ein bestimmtes Pseudopotential erzwingen, z. B. Fe=Fe.rel-pbe.UPF. Kann wiederholt werden. Ohne diese Angabe wählt Olla-DFT mit 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | Cutoff der Wellenfunktionen in Ry (ohne Angabe der von den UPF-Dateien empfohlene) |
| `--ecutrho ECUTRHO` | Cutoff der Dichte in Ry (ohne Angabe der von den UPF-Dateien empfohlene) |
| `--pw-cmd` | pw.x-Programm für --run; aus seinem Pfad werden die übrigen QE-Binärdateien gefunden |
| `--nproc NPROC` | Anzahl der MPI-Prozesse für die mit --run gestarteten Rechnungen |
| `--timeout TIMEOUT` | Zeitlimit in Sekunden für jeden pw.x-Lauf |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria effmass`](THEORY.md)

### `wannier`

Wannier-Funktionen: Bänder, Zentren und Ausdehnung interpolieren, ohne wannier90 zu brauchen

**Verwendung:** `olla-dft wannier [-h] [-o OUTDIR] [-g NxNxN] [-p SITIO:ORBITAL] [--bands BANDS] [--exclude 5-8] [--window MIN:MAX] [--frozen MIN:MAX] [--no-minimize] [--iterations ITERATIONS] [--points POINTS] [--dft-bands DIR] [--no-dft-bands] [--dos N] [--sigma SIGMA] [--run] [--collect] [--pw-cmd PW_CMD] [--pw2wan-cmd PW2WAN_CMD] [--nproc NPROC] [--timeout TIMEOUT] [--pseudo-dir PSEUDO_DIR] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kgrid NxNxN] [--insulator] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] [file]`

**Argumente:**

- `file` — Struktur

**Optionen:**

| Option | Beschreibung |
|---|---|
| `-o, --outdir` | Ausgabeordner (Standard: `wannier`) |
| `-g, --grid NxNxN` | VOLLSTÄNDIGES k-Punkt-Gitter (Standard 4x4x4). Es bestimmt die Qualität der Interpolation |
| `-p, --projections SITIO:ORBITAL` | Testorbitale: 'Si:sp3', 'O:p;Ti:d', 'f=0.125,0.125,0.125:s'. Mehrere durch ';' getrennt. Mit 'auto' werden s und p auf jedes Atom gesetzt (Standard: `auto`) |
| `--bands BANDS` | Bänder des nscf (Standard: so viele wie nötig) |
| `--exclude 5-8` | Bänder, die NICHT in die Wannierisierung eingehen |
| `--window MIN:MAX` | äußeres Entflechtungsfenster in eV: aus welchen Bändern der Unterraum gewählt werden kann. Nötig, wenn die Bänder mit anderen verflochten sind (Leitung, Metalle) |
| `--frozen MIN:MAX` | eingefrorenes Fenster in eV: die Bänder darin werden EXAKT reproduziert. Meist die Valenz plus das Stück Leitungsband, das dich interessiert |
| `--no-minimize` | in der Projektions-Gauge bleiben, ohne die Ausdehnung zu minimieren |
| `--iterations ITERATIONS` | Minimierungsschritte (Standard 500) |
| `--points POINTS` | Punkte pro Abschnitt des interpolierten Pfads (Standard: `30`) |
| `--dft-bands DIR` | Ordner mit der DFT-Bandrechnung zum Vergleich; ohne ihn gibt es keine echte Validierung |
| `--no-dft-bands` | mit --run Schritt 4 (Bänder) überspringen |
| `--dos N` | zusätzlich interpolierte DOS auf einem NxNxN-Gitter |
| `--sigma SIGMA` | Verbreiterung der interpolierten DOS (eV) (Standard: `0.05`) |
| `--run` | die vier Schritte der Reihe nach starten |
| `--collect` | analysieren, was bereits gelaufen ist |
| `--pw-cmd` | pw.x-Programm für --run; aus seinem Pfad werden die übrigen QE-Binärdateien gefunden |
| `--pw2wan-cmd` | pw2wannier90.x-Programm (Standard: neben pw.x) |
| `--nproc NPROC` | Anzahl der MPI-Prozesse für die mit --run gestarteten Rechnungen |
| `--timeout TIMEOUT` | Zeitlimit in Sekunden für jeden pw.x-Lauf |
| `--pseudo-dir` | Ordner mit den UPF-Pseudopotentialen (ohne Angabe der aus 'olla-dft config') |
| `--ecutwfc ECUTWFC` | Cutoff der Wellenfunktionen in Ry (ohne Angabe der von den UPF-Dateien empfohlene) |
| `--ecutrho ECUTRHO` | Cutoff der Dichte in Ry (ohne Angabe der von den UPF-Dateien empfohlene) |
| `--kgrid NxNxN` | Gitter des anfänglichen scf |
| `--insulator` | occupations='fixed' (Isolatoren; Standard: smearing) |
| `--dpi DPI` | Auflösung der Bitmap-Abbildungen in Punkten pro Zoll (Standard: 600) |
| `--format` | kommagetrennte Abbildungsformate: pdf,png,svg,eps,tif (Standard: pdf,png) |
| `--no-plot` | nur die Daten exportieren, ohne die Abbildung zu erzeugen |
| `-t, --template` | visuelle Vorlage der Abbildung (die Liste zeigt 'olla-dft templates list') |
| `--font {sans,serif,latex}` | Schriftfamilie: sans, serif oder latex (Computer Modern) |
| `--usetex` | die Texte der Abbildung mit echtem LaTeX setzen (braucht eine LaTeX-Installation) |
| `--palette` | Farbpalette: grayscale, okabe-ito, okabe-ito-dark oder kommagetrennte Hexfarben |
| `--background` | Hintergrundfarbe der Abbildung, z. B. '#FFFFFF' oder 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | Abbildungsbreite je nach Zeitschrift (Standard: generic) |
| `--width` | Breite der Abbildung: single, onehalf oder double, oder eine Zahl in Millimetern |
| `--mono` | Graustufenversion: schwarze Tinte und Linienmuster |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria wannier`](THEORY.md)

### `unfold`

die Bänder einer Superzelle auf die primitive Brillouin-Zone entfalten

**Verwendung:** `olla-dft unfold [-h] [-o OUTDIR] [--prefix PREFIX] [--bands BANDS] [--spin {up,dw}] [--emin EMIN] [--emax EMAX] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] path primitive`

**Argumente:**

- `path` — Ordner der Bandrechnung der Superzelle
- `primitive` — Struktur der PRIMITIVEN Zelle

**Optionen:**

| Option | Beschreibung |
|---|---|
| `-o, --outdir` | Ausgabeordner (Standard: `.`) |
| `--prefix` | Prefix der Rechnung (wird automatisch erkannt) |
| `--bands BANDS` | wie viele Bänder entfaltet werden (ab dem niedrigsten) |
| `--spin {up,dw}` | zu entfaltender Spinkanal, wenn die Rechnung lsda ist (EIN Kanal wird pro Lauf entfaltet; standardmäßig up) (Standard: `up`) |
| `--emin EMIN` | untere Grenze der Energieachse (eV) (Standard: `-6.0`) |
| `--emax EMAX` | obere Grenze der Energieachse (eV) (Standard: `6.0`) |
| `--dpi DPI` | Auflösung der Bitmap-Abbildungen in Punkten pro Zoll (Standard: 600) |
| `--format` | kommagetrennte Abbildungsformate: pdf,png,svg,eps,tif (Standard: pdf,png) |
| `--no-plot` | nur die Daten exportieren, ohne die Abbildung zu erzeugen |
| `-t, --template` | visuelle Vorlage der Abbildung (die Liste zeigt 'olla-dft templates list') |
| `--font {sans,serif,latex}` | Schriftfamilie: sans, serif oder latex (Computer Modern) |
| `--usetex` | die Texte der Abbildung mit echtem LaTeX setzen (braucht eine LaTeX-Installation) |
| `--palette` | Farbpalette: grayscale, okabe-ito, okabe-ito-dark oder kommagetrennte Hexfarben |
| `--background` | Hintergrundfarbe der Abbildung, z. B. '#FFFFFF' oder 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | Abbildungsbreite je nach Zeitschrift (Standard: generic) |
| `--width` | Breite der Abbildung: single, onehalf oder double, oder eine Zahl in Millimetern |
| `--mono` | Graustufenversion: schwarze Tinte und Linienmuster |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria unfold`](THEORY.md)

### `topology`

Chern-Zahl und Wilson-Schleifen eines Wannier-Modells

**Verwendung:** `olla-dft topology [-h] (--occupied N | --fermi EV) [-g NxN] [--plane {xy,xz,yz}] [--fixed K] [--gap-tol EV] [-o OUTDIR] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] MODELO`

**Argumente:**

- `MODELO` — *_hr.dat-Datei oder Ordner, der WANNIER_hr.dat enthält

**Optionen:**

| Option | Beschreibung |
|---|---|
| `--occupied N` | Anzahl der besetzten Bänder des isolierten Unterraums |
| `--fermi EV` | Fermi-Niveau; wird abgelehnt, wenn es ein Band schneidet |
| `-g, --grid NxN` | periodisches Gitter des 2D-Schnitts (Standard: 40x40) |
| `--plane {xy,xz,yz}` | orientierte Ebene des BZ-Schnitts (Standard: xy) |
| `--fixed K` | senkrechte fraktionelle Koordinate (Standard: 0) |
| `--gap-tol EV` | minimale direkte Bandlücke, um die Invariante zu akzeptieren (Standard: 1e-8) |
| `-o, --outdir` | Ausgabeordner (Standard: `topology`) |
| `--dpi DPI` | Auflösung der Bitmap-Abbildungen in Punkten pro Zoll (Standard: 600) |
| `--format` | kommagetrennte Abbildungsformate: pdf,png,svg,eps,tif (Standard: pdf,png) |
| `--no-plot` | nur die Daten exportieren, ohne die Abbildung zu erzeugen |
| `-t, --template` | visuelle Vorlage der Abbildung (die Liste zeigt 'olla-dft templates list') |
| `--font {sans,serif,latex}` | Schriftfamilie: sans, serif oder latex (Computer Modern) |
| `--usetex` | die Texte der Abbildung mit echtem LaTeX setzen (braucht eine LaTeX-Installation) |
| `--palette` | Farbpalette: grayscale, okabe-ito, okabe-ito-dark oder kommagetrennte Hexfarben |
| `--background` | Hintergrundfarbe der Abbildung, z. B. '#FFFFFF' oder 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | Abbildungsbreite je nach Zeitschrift (Standard: generic) |
| `--width` | Breite der Abbildung: single, onehalf oder double, oder eine Zahl in Millimetern |
| `--mono` | Graustufenversion: schwarze Tinte und Linienmuster |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria topology`](THEORY.md)

### `hubbard`

Hubbard-U per linearer Antwort (hp.x), statt es aus einem Artikel zu kopieren

**Verwendung:** `olla-dft hubbard [-h] [-o OUTDIR] [--species SPECIES] [--qgrid QGRID] [--projection {atomic,ortho-atomic,norm-atomic,wannier,pseudo}] [--hubbard-style {legacy,card}] [--cycle] [--max-iter MAX_ITER] [--tol TOL] [--mixing MIXING] [--collect] [--pw-cmd PW_CMD] [--nproc NPROC] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--metal] [--nspin {1,2}] [--mag MAG] [--intersite] [--v-threshold eV] file`

**Argumente:**

- `file` — Struktur

**Optionen:**

| Option | Beschreibung |
|---|---|
| `-o, --outdir` | Ausgabeordner (Standard: `hubbard`) |
| `--species` | zu störende Spezies, kommagetrennt. Standardmäßig die Übergangsmetalle und Seltenen Erden der Struktur |
| `--qgrid` | q-Gitter der linearen Antwort; entspricht einer Superzelle aus nq1*nq2*nq3 Zellen (Standard: `2x2x2`) |
| `--projection {atomic,ortho-atomic,norm-atomic,wannier,pseudo}` | Projektionsschema. Das U gilt NUR mit demselben Schema, mit dem es berechnet wurde (Standard: `ortho-atomic`) |
| `--hubbard-style {legacy,card}` | DFT+U-Syntax des scf: legacy = lda_plus_u (QE <= 7.0), card = HUBBARD-Karte (QE >= 7.1, wo die alte Syntax ein Fehler ist) (Standard: `legacy`) |
| `--cycle` | vollständiger Selbstkonsistenzzyklus: scf -> hp.x -> scf mit dem neuen U, bis es sich nicht mehr ändert |
| `--max-iter MAX_ITER` | maximale Iterationen des Zyklus scf -> hp.x -> scf mit --cycle (Standard: 6) |
| `--tol TOL` | Änderung in eV, unterhalb derer es als konvergiert gilt (Standard: `0.05`) |
| `--mixing MIXING` | Dämpfung des Schritts; auf 0.5 senken, wenn es oszilliert (Standard: `1.0`) |
| `--collect` | die Ergebnisse einer bereits gelaufenen Rechnung einlesen, statt die Eingaben vorzubereiten |
| `--pw-cmd` | pw.x-Programm für --run; aus seinem Pfad werden die übrigen QE-Binärdateien gefunden |
| `--nproc NPROC` | Anzahl der MPI-Prozesse für die mit --run gestarteten Rechnungen |
| `--pseudo-dir` | Ordner mit den UPF-Pseudopotentialen (ohne Angabe der aus 'olla-dft config') |
| `--pseudo EL=UPF` | ein bestimmtes Pseudopotential erzwingen, z. B. Fe=Fe.rel-pbe.UPF. Kann wiederholt werden. Ohne diese Angabe wählt Olla-DFT mit 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | Cutoff der Wellenfunktionen in Ry (ohne Angabe der von den UPF-Dateien empfohlene) |
| `--ecutrho ECUTRHO` | Cutoff der Dichte in Ry (ohne Angabe der von den UPF-Dateien empfohlene) |
| `--kspacing KSPACING` | Abstand des k-Gitters in Å^-1 |
| `--metal` | metallisches System: Besetzungen mit Smearing statt fester |
| `--nspin {1,2}` | 2 schaltet die Spinpolarisation ein (Standard: 1) |
| `--mag` | Startmagnetisierung: eine Zahl (0.5) oder pro Element (Fe=0.7,O=0). Impliziert --nspin 2 |
| `--intersite` | zusätzlich zu den U-Werten die Intersite-V lesen, die hp.x bereits schreibt, und die HUBBARD-Karte für QE >= 7.1 erzeugen |
| `--v-threshold eV` | V unterhalb dieses Werts wird weder aufgelistet noch geschrieben (Standard: `0.01`) |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria hubbard`](THEORY.md)

## Spektren und Antwort

### `optics`

ε(ω), Absorption und Tauc mit epsilon.x (NC-Pseudopotentiale)

**Verwendung:** `olla-dft optics [-h] [-o OUTDIR] [--run] [--collect] [--pw-cmd PW_CMD] [--nproc NPROC] [-j N] [--redo] [--max-time T] [--estimate] [--timeout TIMEOUT] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--insulator] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--size {paper,poster,presentation}] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--aspect ASPECT] [--mono] [--wmax WMAX] [--smear SMEAR] [--metal] [--suite] [--tauc {direct,indirect}] [--scissor SCISSOR] file`

**Argumente:**

- `file` — Struktur (CIF, POSCAR, pw.x-Eingabe...)

**Optionen:**

| Option | Beschreibung |
|---|---|
| `-o, --outdir` | Ordner des Sweeps (Standard: `opticas`) |
| `--run` | die Rechnungen jetzt ausführen, eine nach der anderen |
| `--collect` | nur bereits gelaufene Rechnungen analysieren |
| `--pw-cmd` | pw.x-Programm (überschreibt die Konfiguration) |
| `--nproc NPROC` | MPI-Prozesse pro Rechnung |
| `-j, --jobs N` | gleichzeitige Rechnungen (Standard: 1). Ohne --nproc werden die Threads der Maschine unter ihnen aufgeteilt |
| `--redo` | auch bereits abgeschlossene Rechnungen wiederholen |
| `--max-time T` | GESAMTES Zeitbudget: 90m, 2h, 3600. Ist es aufgebraucht, werden keine weiteren gestartet und der Sweep bleibt fortsetzbar |
| `--estimate` | abschätzen, wie lange der Sweep dauern wird, und beenden, anhand der Historie von 'olla-dft db' |
| `--timeout TIMEOUT` | Limit in Sekunden pro Rechnung |
| `--pseudo-dir` | Pseudopotential-Ordner |
| `--pseudo EL=UPF` | ein bestimmtes Pseudopotential erzwingen, z. B. Fe=Fe.rel-pbe.UPF. Kann wiederholt werden. Ohne diese Angabe wählt Olla-DFT mit 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | Wellen-Cutoff (Ry) |
| `--ecutrho ECUTRHO` | Cutoff der Dichte (Ry) |
| `--kspacing KSPACING` | k-Abstand in Å^-1 |
| `--insulator` | occupations='fixed' |
| `--dpi DPI` | Auflösung der Bitmap-Abbildungen in Punkten pro Zoll (Standard: 600) |
| `--format` | kommagetrennte Abbildungsformate: pdf,png,svg,eps,tif (Standard: pdf,png) |
| `--no-plot` | nur die Daten exportieren, ohne die Abbildung zu erzeugen |
| `-t, --template` | visuelle Vorlage der Abbildung |
| `--size {paper,poster,presentation}` | Abbildungsgröße: paper, presentation oder poster |
| `--font {sans,serif,latex}` | Schriftfamilie: sans, serif oder latex (Computer Modern) |
| `--usetex` | die Texte der Abbildung mit echtem LaTeX setzen (braucht eine LaTeX-Installation) |
| `--palette` | Farbpalette: grayscale, okabe-ito, okabe-ito-dark oder kommagetrennte Hexfarben |
| `--background` | Hintergrundfarbe der Abbildung, z. B. '#FFFFFF' oder 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | Abbildungsbreite je nach Zeitschrift (Standard: generic) |
| `--width` | Breite der Abbildung: single, onehalf oder double, oder eine Zahl in Millimetern |
| `--aspect ASPECT` | Verhältnis Höhe/Breite der Abbildung |
| `--mono` | Graustufenversion: schwarze Tinte und Linienmuster |
| `--wmax WMAX` | maximale Energie des Spektrums (eV) (Standard: `20.0`) |
| `--smear SMEAR` | Interband-Verbreiterung (eV) (Standard: `0.1`) |
| `--metal` | metallisches System (Besetzungen mit Smearing) |
| `--suite` | zusätzlich ein Austausch-JSON für die anderen Apps der Suite exportieren |
| `--tauc {direct,indirect}` | Übergangstyp für den Tauc-Plot (Standard: `direct`) |
| `--scissor SCISSOR` | starre Verschiebung der Bandlücke in eV (experimentelle oder GW-Lücke minus berechnete Lücke); verschiebt ε2 und baut ε1 per Kramers-Kronig neu auf |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria optics`](THEORY.md)

### `tddft`

optische Absorption mit TDDFPT: lässt das angeregte Elektron und sein Loch einander sehen

**Verwendung:** `olla-dft tddft [-h] [-o OUTDIR] [--method {lanczos,davidson}] [--iter ITER] [--pol {1,2,3,4}] [--states STATES] [--emin EMIN] [--emax EMAX] [--broadening BROADENING] [--scissor SCISSOR] [--extrapolation {no,constant,osc}] [--tamm-dancoff] [--rpa] [--gamma] [--gap GAP] [--compare OPTICS.dat] [--nbnd NBND] [--collect] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--metal] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] [file]`

**Argumente:**

- `file` — Struktur

**Optionen:**

| Option | Beschreibung |
|---|---|
| `-o, --outdir` | Ausgabeordner (Standard: `tddft`) |
| `--method {lanczos,davidson}` | lanczos liefert das ganze Spektrum; davidson liefert die ersten Anregungen eine nach der anderen (Standard: `lanczos`) |
| `--iter ITER` | Lanczos-Iterationen: sie bestimmen die Auflösung (Standard: `500`) |
| `--pol {1,2,3,4}` | 1/2/3 = xx/yy/zz, 4 = voller Tensor (Standard: `4`) |
| `--states STATES` | zu suchende Anregungen (davidson) (Standard: `10`) |
| `--emin EMIN` | untere Grenze der Energieachse (eV) |
| `--emax EMAX` | obere Grenze der Energieachse (eV) (Standard: `15.0`) |
| `--broadening BROADENING` | Verbreiterung in eV (Standard 0.05). Mit --collect legt sie die Nachweisschwelle des Exzitons fest; ohne Angabe wird sie aus spectrum.in gelesen |
| `--scissor SCISSOR` | starre Verschiebung der leeren Bänder in eV (nur lanczos): kompensiert die unterschätzte Bandlücke |
| `--extrapolation {no,constant,osc}` | Extrapolation der Lanczos-Kette im Spektrum: no, constant oder osc (Standard: osc) |
| `--tamm-dancoff` | Tamm-Dancoff-Näherung: billiger, nicht exakt |
| `--rpa` | den xc-Kernel abschalten, um zu sehen, wie viel er beiträgt |
| `--gamma` | K_POINTS gamma erzwingen. Wird automatisch erkannt, wenn die Struktur ein Molekül ist |
| `--gap GAP` | Bandlücke unabhängiger Teilchen in eV, um zu erkennen, ob ein gebundenes Exziton vorliegt |
| `--compare OPTICS.dat` | das Spektrum von 'olla-dft optics' überlagern |
| `--nbnd NBND` | Anzahl der Bänder des scf (standardmäßig, was pw.x entscheidet) |
| `--collect` | die Ergebnisse einer bereits gelaufenen Rechnung einlesen, statt die Eingaben vorzubereiten |
| `--pseudo-dir` | Ordner mit den UPF-Pseudopotentialen (ohne Angabe der aus 'olla-dft config') |
| `--pseudo EL=UPF` | ein bestimmtes Pseudopotential erzwingen, z. B. Fe=Fe.rel-pbe.UPF; kann wiederholt werden |
| `--ecutwfc ECUTWFC` | Cutoff der Wellenfunktionen in Ry (ohne Angabe der von den UPF-Dateien empfohlene) |
| `--ecutrho ECUTRHO` | Cutoff der Dichte in Ry (ohne Angabe der von den UPF-Dateien empfohlene) |
| `--kspacing KSPACING` | Abstand des k-Gitters in Å^-1 |
| `--metal` | metallisches System: Besetzungen mit Smearing statt fester |
| `--dpi DPI` | Auflösung der Bitmap-Abbildungen in Punkten pro Zoll (Standard: 600) |
| `--format` | kommagetrennte Abbildungsformate: pdf,png,svg,eps,tif (Standard: pdf,png) |
| `--no-plot` | nur die Daten exportieren, ohne die Abbildung zu erzeugen |
| `-t, --template` | visuelle Vorlage der Abbildung (die Liste zeigt 'olla-dft templates list') |
| `--font {sans,serif,latex}` | Schriftfamilie: sans, serif oder latex (Computer Modern) |
| `--usetex` | die Texte der Abbildung mit echtem LaTeX setzen (braucht eine LaTeX-Installation) |
| `--palette` | Farbpalette: grayscale, okabe-ito, okabe-ito-dark oder kommagetrennte Hexfarben |
| `--background` | Hintergrundfarbe der Abbildung, z. B. '#FFFFFF' oder 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | Abbildungsbreite je nach Zeitschrift (Standard: generic) |
| `--width` | Breite der Abbildung: single, onehalf oder double, oder eine Zahl in Millimetern |
| `--mono` | Graustufenversion: schwarze Tinte und Linienmuster |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria tddft`](THEORY.md)

### `xanes`

XANES/NEXAFS: Röntgenabsorption nahe der Kante (xspectra.x)

**Verwendung:** `olla-dft xanes [-h] [-o OUTDIR] [--element ELEMENT] [--site SITE] [--edge EDGE] [--core-hole UPF] [--polarization POLARIZATION] [--average] [--emin EMIN] [--emax EMAX] [--broadening BROADENING] [--r-paw R_PAW] [--collect] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--metal] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] file`

**Argumente:**

- `file` — Struktur

**Optionen:**

| Option | Beschreibung |
|---|---|
| `-o, --outdir` | Ausgabeordner (Standard: `xanes`) |
| `--element` | absorbierendes Element |
| `--site SITE` | welches Atom dieses Elements (ab 0) |
| `--edge` | Kante: K, L1, L2, L3 oder L23 (die, die xspectra.x berechnet; M-Kanten nicht) (Standard: `K`) |
| `--core-hole UPF` | Pseudopotential mit Core-Loch (olla-dft corehole) |
| `--polarization` | Richtung des elektrischen Feldes, z. B. '0 0 1' (Standard: `1 0 0`) |
| `--average` | drei orthogonale Richtungen und ihr Mittel: das entspricht einer Pulverprobe |
| `--emin EMIN` | untere Grenze der Energieachse (eV) (Standard: `-10.0`) |
| `--emax EMAX` | obere Grenze der Energieachse (eV) (Standard: `30.0`) |
| `--broadening BROADENING` | Verbreiterung in eV (xgamma) (Standard: `0.8`) |
| `--r-paw R_PAW` | Radius der PAW-Kugel des Absorbers für xspectra.x, in bohr (Standard: 3.0) |
| `--collect` | die Ergebnisse einer bereits gelaufenen Rechnung einlesen, statt die Eingaben vorzubereiten |
| `--pseudo-dir` | Ordner mit den UPF-Pseudopotentialen (ohne Angabe der aus 'olla-dft config') |
| `--pseudo EL=UPF` | ein bestimmtes Pseudopotential erzwingen, z. B. Fe=Fe.rel-pbe.UPF. Kann wiederholt werden. Ohne diese Angabe wählt Olla-DFT mit 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | Cutoff der Wellenfunktionen in Ry (ohne Angabe der von den UPF-Dateien empfohlene) |
| `--ecutrho ECUTRHO` | Cutoff der Dichte in Ry (ohne Angabe der von den UPF-Dateien empfohlene) |
| `--kspacing KSPACING` | Abstand des k-Gitters in Å^-1 |
| `--metal` | metallisches System: Besetzungen mit Smearing statt fester |
| `--dpi DPI` | Auflösung der Bitmap-Abbildungen in Punkten pro Zoll (Standard: 600) |
| `--format` | kommagetrennte Abbildungsformate: pdf,png,svg,eps,tif (Standard: pdf,png) |
| `--no-plot` | nur die Daten exportieren, ohne die Abbildung zu erzeugen |
| `-t, --template` | visuelle Vorlage der Abbildung (die Liste zeigt 'olla-dft templates list') |
| `--font {sans,serif,latex}` | Schriftfamilie: sans, serif oder latex (Computer Modern) |
| `--usetex` | die Texte der Abbildung mit echtem LaTeX setzen (braucht eine LaTeX-Installation) |
| `--palette` | Farbpalette: grayscale, okabe-ito, okabe-ito-dark oder kommagetrennte Hexfarben |
| `--background` | Hintergrundfarbe der Abbildung, z. B. '#FFFFFF' oder 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | Abbildungsbreite je nach Zeitschrift (Standard: generic) |
| `--width` | Breite der Abbildung: single, onehalf oder double, oder eine Zahl in Millimetern |
| `--mono` | Graustufenversion: schwarze Tinte und Linienmuster |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria xanes`](THEORY.md)

### `xps`

Core-Level-Verschiebungen (Anfangszustand)

**Verwendung:** `olla-dft xps [-h] [-o OUTDIR] [--core-hole EL=UPF] [--collect] [--suite] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--metal] file`

**Argumente:**

- `file` — Eingabestruktur (CIF, POSCAR, pw.x-Eingabe...)

**Optionen:**

| Option | Beschreibung |
|---|---|
| `-o, --outdir` | Ausgabeordner (Standard: `xps`) |
| `--core-hole EL=UPF` | Pseudopotential mit Core-Loch, z. B. Si=Si.star1s.UPF. Kann wiederholt werden. Ohne dies liefert initial_state.x eine Tabelle voller Nullen |
| `--collect` | die Ergebnisse einer bereits gelaufenen Rechnung einlesen, statt die Eingaben vorzubereiten |
| `--suite` | zusätzlich ein Austausch-JSON für die anderen Apps der Suite exportieren |
| `--pseudo-dir` | Ordner mit den UPF-Pseudopotentialen (ohne Angabe der aus 'olla-dft config') |
| `--pseudo EL=UPF` | ein bestimmtes Pseudopotential erzwingen, z. B. Fe=Fe.rel-pbe.UPF. Kann wiederholt werden. Ohne diese Angabe wählt Olla-DFT mit 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | Cutoff der Wellenfunktionen in Ry (ohne Angabe der von den UPF-Dateien empfohlene) |
| `--ecutrho ECUTRHO` | Cutoff der Dichte in Ry (ohne Angabe der von den UPF-Dateien empfohlene) |
| `--kspacing KSPACING` | Abstand des k-Gitters in Å^-1 |
| `--metal` | metallisches System: Besetzungen mit Smearing statt fester |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria xps`](THEORY.md)

### `corehole`

das Pseudopotentialpaar normal + Core-Loch erzeugen (ld1.x), für XPS und XANES

**Verwendung:** `olla-dft corehole [-h] [-o OUTDIR] [--edge EDGE] [--functional FUNCTIONAL] [--rcut RCUT] [--rel {0,1,2}] [--semicore] [--pseudotype {1,2,3}] [--plain] [--only-inputs] [--projectors {1,2}] [--ld1-cmd LD1_CMD] [--core-wfc UPF] [--orbital ORBITAL] [--output OUTPUT] [element]`

**Argumente:**

- `element` — Elementsymbol, z. B. Si

**Optionen:**

| Option | Beschreibung |
|---|---|
| `-o, --outdir` | Ausgabeordner (Standard: `pseudos`) |
| `--edge` | Kante/Niveau des Lochs: K (1s), L1 (2s), L23 (2p), M1, M23, M45 (Standard: `K`) |
| `--functional` | Funktional des Pseudopotentials; es muss dasselbe sein, mit dem du pw.x ausführen wirst (Standard: `PBE`) |
| `--rcut RCUT` | Cutoff-Radius in bohr (standardmäßig einer pro Zeile des Periodensystems) |
| `--rel {0,1,2}` | 0 nichtrelativistisch, 1 skalar, 2 voll |
| `--semicore` | die Schale (n-1)s(n-1)p in die Valenz aufnehmen |
| `--pseudotype {1,2,3}` | 1 und 2 sind normerhaltend, 3 ist ultrasoft (Standard: `2`) |
| `--plain` | NUR das normale Pseudopotential erzeugen, ohne das mit Core-Loch. Nützlich, um ein konsistentes Pseudo für ein Element zu bekommen, das es nicht unterstützt |
| `--only-inputs` | die ld1.x-Eingaben schreiben, ohne sie auszuführen |
| `--projectors {1,2}` | GIPAW-Projektoren pro Kanal. XSpectra empfiehlt 2, aber mit 2 wird das Pseudo ultrasoft und --rcut muss meist von Hand angepasst werden (Standard: `1`) |
| `--ld1-cmd` | Pfad zu ld1.x |
| `--core-wfc UPF` | statt zu erzeugen: die Core-Wellenfunktion aus einer UPF-Datei in dem Format extrahieren, das xspectra.x liest |
| `--orbital` | zu prüfendes Orbital, z. B. 1S |
| `--output` | Ausgabedatei für --core-wfc |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria corehole`](THEORY.md)

### `charge`

Ladungsdichte / ELF / Spin mit pp.x

**Verwendung:** `olla-dft charge [-h] [-o OUTDIR] [--field {density,elf,spin,potential,vtotal}] [--axis AXIS] [--rerun] [--pw-cmd PW_CMD] [--nproc NPROC] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] [path]`

**Argumente:**

- `path` — Ordner der Rechnung

**Optionen:**

| Option | Beschreibung |
|---|---|
| `-o, --outdir` | Ausgabeordner (Standard: `.`) |
| `--field {density,elf,spin,potential,vtotal}` | mit pp.x zu berechnendes Feld: density, elf, spin, potential oder vtotal (Standard: density) |
| `--axis` | Achse des planaren Profils (Standard: `c`) |
| `--rerun` | pp.x erneut ausführen, auch wenn die Cube-Datei bereits vorhanden ist |
| `--pw-cmd` | pw.x-Programm für --run; aus seinem Pfad werden die übrigen QE-Binärdateien gefunden |
| `--nproc NPROC` | Anzahl der MPI-Prozesse für die mit --run gestarteten Rechnungen |
| `--dpi DPI` | Auflösung der Bitmap-Abbildungen in Punkten pro Zoll (Standard: 600) |
| `--format` | kommagetrennte Abbildungsformate: pdf,png,svg,eps,tif (Standard: pdf,png) |
| `--no-plot` | nur die Daten exportieren, ohne die Abbildung zu erzeugen |
| `-t, --template` | visuelle Vorlage der Abbildung (die Liste zeigt 'olla-dft templates list') |
| `--font {sans,serif,latex}` | Schriftfamilie: sans, serif oder latex (Computer Modern) |
| `--usetex` | die Texte der Abbildung mit echtem LaTeX setzen (braucht eine LaTeX-Installation) |
| `--palette` | Farbpalette: grayscale, okabe-ito, okabe-ito-dark oder kommagetrennte Hexfarben |
| `--background` | Hintergrundfarbe der Abbildung, z. B. '#FFFFFF' oder 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | Abbildungsbreite je nach Zeitschrift (Standard: generic) |
| `--width` | Breite der Abbildung: single, onehalf oder double, oder eine Zahl in Millimetern |
| `--mono` | Graustufenversion: schwarze Tinte und Linienmuster |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria charge`](THEORY.md)

### `charges`

Löwdin/Bader-Ladungen und Dichtedifferenz

**Verwendung:** `olla-dft charges [-h] [--lowdin LOWDIN] [--bader BADER] [--difference CUBE [CUBE ...]] [--pseudo-dir PSEUDO_DIR] [--axis {0,1,2}] [-o OUTDIR] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] [file]`

**Argumente:**

- `file` — Struktur (für Bader)

**Optionen:**

| Option | Beschreibung |
|---|---|
| `--lowdin` | projwfc.x-Ausgabe |
| `--bader` | Dichte-Cube (plot_num=0) |
| `--difference CUBE` | total.cube teil1.cube teil2.cube ... |
| `--pseudo-dir` | Ordner mit den UPF-Dateien der Rechnung: daraus kommt Z_valence für die Spalte 'neta' (netto) (überschreibt config) |
| `--axis {0,1,2}` | Achse des planaren Profils der Dichtedifferenz: 0, 1 oder 2 (Standard: 2) |
| `-o, --outdir` | Ausgabeordner (Standard: `.`) |
| `--dpi DPI` | Auflösung der Bitmap-Abbildungen in Punkten pro Zoll (Standard: 600) |
| `--format` | kommagetrennte Abbildungsformate: pdf,png,svg,eps,tif (Standard: pdf,png) |
| `--no-plot` | nur die Daten exportieren, ohne die Abbildung zu erzeugen |
| `-t, --template` | visuelle Vorlage der Abbildung (die Liste zeigt 'olla-dft templates list') |
| `--font {sans,serif,latex}` | Schriftfamilie: sans, serif oder latex (Computer Modern) |
| `--usetex` | die Texte der Abbildung mit echtem LaTeX setzen (braucht eine LaTeX-Installation) |
| `--palette` | Farbpalette: grayscale, okabe-ito, okabe-ito-dark oder kommagetrennte Hexfarben |
| `--background` | Hintergrundfarbe der Abbildung, z. B. '#FFFFFF' oder 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | Abbildungsbreite je nach Zeitschrift (Standard: generic) |
| `--width` | Breite der Abbildung: single, onehalf oder double, oder eine Zahl in Millimetern |
| `--mono` | Graustufenversion: schwarze Tinte und Linienmuster |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria charges`](THEORY.md)

### `wf`

Austrittsarbeit aus einer Rechnung mit Vakuum

**Verwendung:** `olla-dft wf [-h] [-o OUTDIR] [--axis AXIS] [--rerun] [--pw-cmd PW_CMD] [--nproc NPROC] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] [path]`

**Argumente:**

- `path` — Ordner der Rechnung

**Optionen:**

| Option | Beschreibung |
|---|---|
| `-o, --outdir` | Ausgabeordner (Standard: `.`) |
| `--axis` | Vakuumachse: a/b/c (Standard c) |
| `--rerun` | pp.x erneut ausführen, auch wenn die Cube-Datei existiert |
| `--pw-cmd` | pw.x-Programm für --run; aus seinem Pfad werden die übrigen QE-Binärdateien gefunden |
| `--nproc NPROC` | Anzahl der MPI-Prozesse für die mit --run gestarteten Rechnungen |
| `--dpi DPI` | Auflösung der Bitmap-Abbildungen in Punkten pro Zoll (Standard: 600) |
| `--format` | kommagetrennte Abbildungsformate: pdf,png,svg,eps,tif (Standard: pdf,png) |
| `--no-plot` | nur die Daten exportieren, ohne die Abbildung zu erzeugen |
| `-t, --template` | visuelle Vorlage der Abbildung (die Liste zeigt 'olla-dft templates list') |
| `--font {sans,serif,latex}` | Schriftfamilie: sans, serif oder latex (Computer Modern) |
| `--usetex` | die Texte der Abbildung mit echtem LaTeX setzen (braucht eine LaTeX-Installation) |
| `--palette` | Farbpalette: grayscale, okabe-ito, okabe-ito-dark oder kommagetrennte Hexfarben |
| `--background` | Hintergrundfarbe der Abbildung, z. B. '#FFFFFF' oder 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | Abbildungsbreite je nach Zeitschrift (Standard: generic) |
| `--width` | Breite der Abbildung: single, onehalf oder double, oder eine Zahl in Millimetern |
| `--mono` | Graustufenversion: schwarze Tinte und Linienmuster |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria wf`](THEORY.md)

### `berry`

Polarisation per Berry-Phase: spontanes ΔP, Born-Ladungen und Ferroelektrizität

**Verwendung:** `olla-dft berry [-h] [-o OUTDIR] [--gdir {1,2,3}] [--nppstr NPPSTR] [--kperp NxN] [-r ARCHIVO] [--displace ATOMO:dx,dy,dz] [--nlambda NLAMBDA] [--run] [--collect] [--redo] [--pw-cmd PW_CMD] [--nproc NPROC] [--timeout TIMEOUT] [--pseudo-dir PSEUDO_DIR] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kgrid NxNxN] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] file`

**Argumente:**

- `file` — Struktur (die polare, falls es einen Pfad gibt)

**Optionen:**

| Option | Beschreibung |
|---|---|
| `-o, --outdir` | Ausgabeordner (Standard: `berry`) |
| `--gdir {1,2,3}` | Richtung: Vektor des reziproken Gitters (Standard 3) |
| `--nppstr NPPSTR` | Punkte pro k-String (Standard 9); erhöhen, bis sich die Phase nicht mehr ändert |
| `--kperp NxN` | Gitter senkrecht zum String (Standard 6x6) |
| `-r, --reference ARCHIVO` | Referenzstruktur, normalerweise die zentrosymmetrische: ein adiabatischer Pfad zur polaren wird interpoliert und ΔP ist die spontane Polarisation |
| `--displace ATOMO:dx,dy,dz` | Verschiebungspfad eines Atoms, in Å; die Steigung von P liefert die effektive Born-Ladung |
| `--nlambda NLAMBDA` | Punkte des Pfads (Standard 5) |
| `--run` | die Rechnungen ausführen, sobald die Eingaben vorbereitet sind |
| `--collect` | die Ergebnisse einer bereits gelaufenen Rechnung einlesen, statt die Eingaben vorzubereiten |
| `--redo` | auch bereits abgeschlossene Rechnungen wiederholen |
| `--pw-cmd` | pw.x-Programm für --run; aus seinem Pfad werden die übrigen QE-Binärdateien gefunden |
| `--nproc NPROC` | Anzahl der MPI-Prozesse für die mit --run gestarteten Rechnungen |
| `--timeout TIMEOUT` | Zeitlimit in Sekunden für jeden pw.x-Lauf |
| `--pseudo-dir` | Ordner mit den UPF-Pseudopotentialen (ohne Angabe der aus 'olla-dft config') |
| `--ecutwfc ECUTWFC` | Cutoff der Wellenfunktionen in Ry (ohne Angabe der von den UPF-Dateien empfohlene) |
| `--ecutrho ECUTRHO` | Cutoff der Dichte in Ry (ohne Angabe der von den UPF-Dateien empfohlene) |
| `--kgrid NxNxN` | k-Gitter des scf, z. B. 6x6x6 (standardmäßig aus dem k-Abstand) |
| `--dpi DPI` | Auflösung der Bitmap-Abbildungen in Punkten pro Zoll (Standard: 600) |
| `--format` | kommagetrennte Abbildungsformate: pdf,png,svg,eps,tif (Standard: pdf,png) |
| `--no-plot` | nur die Daten exportieren, ohne die Abbildung zu erzeugen |
| `-t, --template` | visuelle Vorlage der Abbildung (die Liste zeigt 'olla-dft templates list') |
| `--font {sans,serif,latex}` | Schriftfamilie: sans, serif oder latex (Computer Modern) |
| `--usetex` | die Texte der Abbildung mit echtem LaTeX setzen (braucht eine LaTeX-Installation) |
| `--palette` | Farbpalette: grayscale, okabe-ito, okabe-ito-dark oder kommagetrennte Hexfarben |
| `--background` | Hintergrundfarbe der Abbildung, z. B. '#FFFFFF' oder 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | Abbildungsbreite je nach Zeitschrift (Standard: generic) |
| `--width` | Breite der Abbildung: single, onehalf oder double, oder eine Zahl in Millimetern |
| `--mono` | Graustufenversion: schwarze Tinte und Linienmuster |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria berry`](THEORY.md)

## Phononen, Transport und Temperatur

### `phonons`

DFPT-Phononen: Dispersion, DOS, Thermodynamik, IR

**Verwendung:** `olla-dft phonons [-h] [-o OUTDIR] [--run] [--collect] [--pw-cmd PW_CMD] [--nproc NPROC] [-j N] [--redo] [--max-time T] [--estimate] [--timeout TIMEOUT] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--insulator] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--size {paper,poster,presentation}] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--aspect ASPECT] [--mono] [--qgrid QGRID] [--gamma] [--raman] [--laser LASER] [--suite] [--tscan T1,T2,...] file`

**Argumente:**

- `file` — Struktur (CIF, POSCAR, pw.x-Eingabe...)

**Optionen:**

| Option | Beschreibung |
|---|---|
| `-o, --outdir` | Ordner des Sweeps (Standard: `fonones`) |
| `--run` | die Rechnungen jetzt ausführen, eine nach der anderen |
| `--collect` | nur bereits gelaufene Rechnungen analysieren |
| `--pw-cmd` | pw.x-Programm (überschreibt die Konfiguration) |
| `--nproc NPROC` | MPI-Prozesse pro Rechnung |
| `-j, --jobs N` | gleichzeitige Rechnungen (Standard: 1). Ohne --nproc werden die Threads der Maschine unter ihnen aufgeteilt |
| `--redo` | auch bereits abgeschlossene Rechnungen wiederholen |
| `--max-time T` | GESAMTES Zeitbudget: 90m, 2h, 3600. Ist es aufgebraucht, werden keine weiteren gestartet und der Sweep bleibt fortsetzbar |
| `--estimate` | abschätzen, wie lange der Sweep dauern wird, und beenden, anhand der Historie von 'olla-dft db' |
| `--timeout TIMEOUT` | Limit in Sekunden pro Rechnung |
| `--pseudo-dir` | Pseudopotential-Ordner |
| `--pseudo EL=UPF` | ein bestimmtes Pseudopotential erzwingen, z. B. Fe=Fe.rel-pbe.UPF. Kann wiederholt werden. Ohne diese Angabe wählt Olla-DFT mit 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | Wellen-Cutoff (Ry) |
| `--ecutrho ECUTRHO` | Cutoff der Dichte (Ry) |
| `--kspacing KSPACING` | k-Abstand in Å^-1 |
| `--insulator` | occupations='fixed' |
| `--dpi DPI` | Auflösung der Bitmap-Abbildungen in Punkten pro Zoll (Standard: 600) |
| `--format` | kommagetrennte Abbildungsformate: pdf,png,svg,eps,tif (Standard: pdf,png) |
| `--no-plot` | nur die Daten exportieren, ohne die Abbildung zu erzeugen |
| `-t, --template` | visuelle Vorlage der Abbildung |
| `--size {paper,poster,presentation}` | Abbildungsgröße: paper, presentation oder poster |
| `--font {sans,serif,latex}` | Schriftfamilie: sans, serif oder latex (Computer Modern) |
| `--usetex` | die Texte der Abbildung mit echtem LaTeX setzen (braucht eine LaTeX-Installation) |
| `--palette` | Farbpalette: grayscale, okabe-ito, okabe-ito-dark oder kommagetrennte Hexfarben |
| `--background` | Hintergrundfarbe der Abbildung, z. B. '#FFFFFF' oder 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | Abbildungsbreite je nach Zeitschrift (Standard: generic) |
| `--width` | Breite der Abbildung: single, onehalf oder double, oder eine Zahl in Millimetern |
| `--aspect ASPECT` | Verhältnis Höhe/Breite der Abbildung |
| `--mono` | Graustufenversion: schwarze Tinte und Linienmuster |
| `--qgrid` | q-Gitter, z. B. 2x2x2 |
| `--gamma` | nur Γ mit dynmat.x: Frequenzen und IR-Aktivitäten |
| `--raman` | zusätzlich Raman-Tensoren und -Intensitäten in Γ (lraman; nur normerhaltende Pseudopotentiale, und deutlich teurer) |
| `--laser LASER` | Laserwellenlänge in nm zur Simulation des Raman-Spektrums (Standard: `532.0`) |
| `--suite` | zusätzlich ein Austausch-JSON exportieren (nur mit --gamma) für die FTIR- und Raman-Apps |
| `--tscan T1,T2,...` | Sweep der ELEKTRONISCHEN Temperatur in K: wiederholt die Phononen mit fermi-dirac-Smearing bei jeder und prüft, ob sich eine imaginäre Mode beim Erwärmen stabilisiert (Ladungsdichtewellen, strukturelle Übergänge) |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria phonons`](THEORY.md)

### `elph`

Elektron-Phonon-Kopplung: lambda, Tc und ein echtes tau für den Transport

**Verwendung:** `olla-dft elph [-h] [-o OUTDIR] [--qgrid QGRID] [--kgrid KGRID] [--kgrid-nscf KGRID_NSCF] [--nsigma NSIGMA] [--sigma SIGMA] [--degauss DEGAUSS] [--debye DEBYE] [--collect] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] [file]`

**Argumente:**

- `file` — Struktur

**Optionen:**

| Option | Beschreibung |
|---|---|
| `-o, --outdir` | Ausgabeordner (Standard: `elph`) |
| `--qgrid` | q-Gitter der DFPT, z. B. 2x2x2 (Standard: 2x2x2) |
| `--kgrid` | k-Gitter des scf |
| `--kgrid-nscf` | k-Gitter des dichten nscf; standardmäßig das Doppelte des scf-Gitters, gerundet auf ein Vielfaches des q-Gitters |
| `--nsigma NSIGMA` | wie viele Verbreiterungen ph.x für lambda durchläuft (el_ph_nsigma; Standard: 10) |
| `--sigma SIGMA` | Schrittweite des Verbreiterungs-Sweeps, in Ry (Standard: `0.005`) |
| `--degauss DEGAUSS` | Smearing des scf in Ry (Standard: 0.02) |
| `--debye DEBYE` | Debye-Temperatur in K, um den Bereich zu markieren, in dem die tau-Formel gilt |
| `--collect` | die Ergebnisse einer bereits gelaufenen Rechnung einlesen, statt die Eingaben vorzubereiten |
| `--pseudo-dir` | Ordner mit den UPF-Pseudopotentialen (ohne Angabe der aus 'olla-dft config') |
| `--pseudo EL=UPF` | ein bestimmtes Pseudopotential erzwingen, z. B. Fe=Fe.rel-pbe.UPF. Kann wiederholt werden. Ohne diese Angabe wählt Olla-DFT mit 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | Cutoff der Wellenfunktionen in Ry (ohne Angabe der von den UPF-Dateien empfohlene) |
| `--ecutrho ECUTRHO` | Cutoff der Dichte in Ry (ohne Angabe der von den UPF-Dateien empfohlene) |
| `--dpi DPI` | Auflösung der Bitmap-Abbildungen in Punkten pro Zoll (Standard: 600) |
| `--format` | kommagetrennte Abbildungsformate: pdf,png,svg,eps,tif (Standard: pdf,png) |
| `--no-plot` | nur die Daten exportieren, ohne die Abbildung zu erzeugen |
| `-t, --template` | visuelle Vorlage der Abbildung (die Liste zeigt 'olla-dft templates list') |
| `--font {sans,serif,latex}` | Schriftfamilie: sans, serif oder latex (Computer Modern) |
| `--usetex` | die Texte der Abbildung mit echtem LaTeX setzen (braucht eine LaTeX-Installation) |
| `--palette` | Farbpalette: grayscale, okabe-ito, okabe-ito-dark oder kommagetrennte Hexfarben |
| `--background` | Hintergrundfarbe der Abbildung, z. B. '#FFFFFF' oder 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | Abbildungsbreite je nach Zeitschrift (Standard: generic) |
| `--width` | Breite der Abbildung: single, onehalf oder double, oder eine Zahl in Millimetern |
| `--mono` | Graustufenversion: schwarze Tinte und Linienmuster |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria elph`](THEORY.md)

### `transport`

Seebeck, sigma/tau und Leistungsfaktor (CRTA)

**Verwendung:** `olla-dft transport [-h] [-o OUTDIR] [--run] [--collect] [--pw-cmd PW_CMD] [--nproc NPROC] [-j N] [--redo] [--max-time T] [--estimate] [--timeout TIMEOUT] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--insulator] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--size {paper,poster,presentation}] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--aspect ASPECT] [--mono] [--grid GRID] [--temperatures TEMPERATURES] [--mu-span MU_SPAN] [--metal] [--nspin {1,2}] [--mag MAG] [--spin-resolved] file`

**Argumente:**

- `file` — Struktur (CIF, POSCAR, pw.x-Eingabe...)

**Optionen:**

| Option | Beschreibung |
|---|---|
| `-o, --outdir` | Ordner des Sweeps (Standard: `transporte`) |
| `--run` | die Rechnungen jetzt ausführen, eine nach der anderen |
| `--collect` | nur bereits gelaufene Rechnungen analysieren |
| `--pw-cmd` | pw.x-Programm (überschreibt die Konfiguration) |
| `--nproc NPROC` | MPI-Prozesse pro Rechnung |
| `-j, --jobs N` | gleichzeitige Rechnungen (Standard: 1). Ohne --nproc werden die Threads der Maschine unter ihnen aufgeteilt |
| `--redo` | auch bereits abgeschlossene Rechnungen wiederholen |
| `--max-time T` | GESAMTES Zeitbudget: 90m, 2h, 3600. Ist es aufgebraucht, werden keine weiteren gestartet und der Sweep bleibt fortsetzbar |
| `--estimate` | abschätzen, wie lange der Sweep dauern wird, und beenden, anhand der Historie von 'olla-dft db' |
| `--timeout TIMEOUT` | Limit in Sekunden pro Rechnung |
| `--pseudo-dir` | Pseudopotential-Ordner |
| `--pseudo EL=UPF` | ein bestimmtes Pseudopotential erzwingen, z. B. Fe=Fe.rel-pbe.UPF. Kann wiederholt werden. Ohne diese Angabe wählt Olla-DFT mit 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | Wellen-Cutoff (Ry) |
| `--ecutrho ECUTRHO` | Cutoff der Dichte (Ry) |
| `--kspacing KSPACING` | k-Abstand in Å^-1 |
| `--insulator` | occupations='fixed' |
| `--dpi DPI` | Auflösung der Bitmap-Abbildungen in Punkten pro Zoll (Standard: 600) |
| `--format` | kommagetrennte Abbildungsformate: pdf,png,svg,eps,tif (Standard: pdf,png) |
| `--no-plot` | nur die Daten exportieren, ohne die Abbildung zu erzeugen |
| `-t, --template` | visuelle Vorlage der Abbildung |
| `--size {paper,poster,presentation}` | Abbildungsgröße: paper, presentation oder poster |
| `--font {sans,serif,latex}` | Schriftfamilie: sans, serif oder latex (Computer Modern) |
| `--usetex` | die Texte der Abbildung mit echtem LaTeX setzen (braucht eine LaTeX-Installation) |
| `--palette` | Farbpalette: grayscale, okabe-ito, okabe-ito-dark oder kommagetrennte Hexfarben |
| `--background` | Hintergrundfarbe der Abbildung, z. B. '#FFFFFF' oder 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | Abbildungsbreite je nach Zeitschrift (Standard: generic) |
| `--width` | Breite der Abbildung: single, onehalf oder double, oder eine Zahl in Millimetern |
| `--aspect ASPECT` | Verhältnis Höhe/Breite der Abbildung |
| `--mono` | Graustufenversion: schwarze Tinte und Linienmuster |
| `--grid` | Gitter des nscf, z. B. 16x16x16 |
| `--temperatures` | kommagetrennte Temperaturen in K (Standard: `300`) |
| `--mu-span MU_SPAN` | Bereich des chemischen Potentials um E_F (eV) (Standard: `1.0`) |
| `--metal` | metallisches System: Besetzungen mit Smearing statt fester |
| `--nspin {1,2}` | 2 schaltet die Spinpolarisation in scf und nscf ein (nötig für --spin-resolved) (Standard: `1`) |
| `--mag` | Startmagnetisierung, z. B. Fe=0.7 (impliziert --nspin 2) |
| `--spin-resolved` | die beiden Spinkanäle trennen (Zweistrommodell) und die Polarisation der Leitfähigkeit sowie die Spin-Thermokraft angeben |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria transport`](THEORY.md)

### `ballistic`

ballistischer Landauer-Leitwert (pwcond.x), für Nanokontakte und Moleküle zwischen Elektroden

**Verwendung:** `olla-dft ballistic [-h] [--scatterer SCATTERER] [-o OUTDIR] [--ikind {0,1}] [--emin EMIN] [--emax EMAX] [--points POINTS] [--nz1 NZ1] [--collect] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] [file]`

**Argumente:**

- `file` — Elektrode: die in z periodische Zelle

**Optionen:**

| Option | Beschreibung |
|---|---|
| `--scatterer` | Streuregion (das Molekül oder der Defekt). Ohne sie kommen nur die komplexen Bänder heraus |
| `-o, --outdir` | Ausgabeordner (Standard: `balistico`) |
| `--ikind {0,1}` | 0 = nur komplexe Bänder, 1 = Leitwert mit derselben Elektrode auf beiden Seiten (Standard: 1 mit --scatterer, sonst 0). Unterschiedliche Elektroden (ikind=2 von pwcond.x) werden nicht unterstützt |
| `--emin EMIN` | untere Grenze der Energieachse (eV) (Standard: `-3.0`) |
| `--emax EMAX` | obere Grenze der Energieachse (eV) (Standard: `3.0`) |
| `--points POINTS` | Anzahl der Energien des Leitwert-Sweeps (Standard: 61) |
| `--nz1 NZ1` | z-Unterteilungen jeder pwcond.x-Scheibe (nz1; Standard: 3) |
| `--collect` | die Ergebnisse einer bereits gelaufenen Rechnung einlesen, statt die Eingaben vorzubereiten |
| `--pseudo-dir` | Ordner mit den UPF-Pseudopotentialen (ohne Angabe der aus 'olla-dft config') |
| `--pseudo EL=UPF` | ein bestimmtes Pseudopotential erzwingen, z. B. Fe=Fe.rel-pbe.UPF; kann wiederholt werden |
| `--ecutwfc ECUTWFC` | Cutoff der Wellenfunktionen in Ry (ohne Angabe der von den UPF-Dateien empfohlene) |
| `--ecutrho ECUTRHO` | Cutoff der Dichte in Ry (ohne Angabe der von den UPF-Dateien empfohlene) |
| `--kspacing KSPACING` | Abstand des k-Gitters in Å^-1 |
| `--dpi DPI` | Auflösung der Bitmap-Abbildungen in Punkten pro Zoll (Standard: 600) |
| `--format` | kommagetrennte Abbildungsformate: pdf,png,svg,eps,tif (Standard: pdf,png) |
| `--no-plot` | nur die Daten exportieren, ohne die Abbildung zu erzeugen |
| `-t, --template` | visuelle Vorlage der Abbildung (die Liste zeigt 'olla-dft templates list') |
| `--font {sans,serif,latex}` | Schriftfamilie: sans, serif oder latex (Computer Modern) |
| `--usetex` | die Texte der Abbildung mit echtem LaTeX setzen (braucht eine LaTeX-Installation) |
| `--palette` | Farbpalette: grayscale, okabe-ito, okabe-ito-dark oder kommagetrennte Hexfarben |
| `--background` | Hintergrundfarbe der Abbildung, z. B. '#FFFFFF' oder 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | Abbildungsbreite je nach Zeitschrift (Standard: generic) |
| `--width` | Breite der Abbildung: single, onehalf oder double, oder eine Zahl in Millimetern |
| `--mono` | Graustufenversion: schwarze Tinte und Linienmuster |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria ballistic`](THEORY.md)

### `kappa`

Gitterwärmeleitfähigkeit: fc3, Phonon-Boltzmann-Gleichung und mittlere freie Weglänge

**Verwendung:** `olla-dft kappa [-h] [-o OUTDIR] [--dim NxNxN] [--dim-fc2 NxNxN] [--distance DISTANCE] [--mesh MESH] [--temps TEMPS] [--isotopes] [--grain UM] [--model MODEL] [--collect] [--force] [--metal] [--pseudo-dir PSEUDO_DIR] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] file`

**Argumente:**

- `file` — Struktur (primitive Zelle)

**Optionen:**

| Option | Beschreibung |
|---|---|
| `-o, --outdir` | Ausgabeordner (Standard: `kappa`) |
| `--dim NxNxN` | Superzelle der fc3 (Standard 2x2x2). Das ist der teure Teil: die Zahl der Konfigurationen wächst schnell |
| `--dim-fc2 NxNxN` | GRÖSSERE Superzelle nur für den harmonischen Teil, der billig ist und mehr Reichweite braucht |
| `--distance DISTANCE` | finite Auslenkung in Å (Standard 0.03) |
| `--mesh MESH` | q-Gitter für die Boltzmann-Gleichung (Standard 13) |
| `--temps` | Temperaturen: 100:800:8 oder 300,500,700 (Standard: `100:800:8`) |
| `--isotopes` | Isotopenstreuung mit natürlichen Häufigkeiten hinzufügen (in Si sind es ~10 %%) |
| `--grain UM` | Korngröße in µm: fügt Streuung an Korngrenzen hinzu |
| `--model` | die Kräfte mit einem gelernten Potential (mace, chgnet, m3gnet) statt mit pw.x berechnen: Sekunden statt Stunden, aber der Absolutwert kann weit danebenliegen |
| `--collect` | die bereits berechneten Kräfte einlesen und lösen |
| `--force` | die Eingaben schreiben, auch wenn es sehr viele sind |
| `--metal` | metallisches System (Besetzungen mit Smearing in den scf-Läufen der fc2/fc3). Ohne dies wird occupations='fixed' verwendet, was für Isolatoren richtig ist |
| `--pseudo-dir` | Ordner mit den UPF-Pseudopotentialen (ohne Angabe der aus 'olla-dft config') |
| `--ecutwfc ECUTWFC` | Cutoff der Wellenfunktionen in Ry (ohne Angabe der von den UPF-Dateien empfohlene) |
| `--ecutrho ECUTRHO` | Cutoff der Dichte in Ry (ohne Angabe der von den UPF-Dateien empfohlene) |
| `--kspacing KSPACING` | Abstand des k-Gitters in Å^-1 (Standard: `0.35`) |
| `--dpi DPI` | Auflösung der Bitmap-Abbildungen in Punkten pro Zoll (Standard: 600) |
| `--format` | kommagetrennte Abbildungsformate: pdf,png,svg,eps,tif (Standard: pdf,png) |
| `--no-plot` | nur die Daten exportieren, ohne die Abbildung zu erzeugen |
| `-t, --template` | visuelle Vorlage der Abbildung (die Liste zeigt 'olla-dft templates list') |
| `--font {sans,serif,latex}` | Schriftfamilie: sans, serif oder latex (Computer Modern) |
| `--usetex` | die Texte der Abbildung mit echtem LaTeX setzen (braucht eine LaTeX-Installation) |
| `--palette` | Farbpalette: grayscale, okabe-ito, okabe-ito-dark oder kommagetrennte Hexfarben |
| `--background` | Hintergrundfarbe der Abbildung, z. B. '#FFFFFF' oder 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | Abbildungsbreite je nach Zeitschrift (Standard: generic) |
| `--width` | Breite der Abbildung: single, onehalf oder double, oder eine Zahl in Millimetern |
| `--mono` | Graustufenversion: schwarze Tinte und Linienmuster |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria kappa`](THEORY.md)

### `qha`

quasiharmonisch: thermische Ausdehnung und a(T)

**Verwendung:** `olla-dft qha [-h] [-o OUTDIR] [--natoms NATOMS] [--cells CELLS] [--cubic] [--structure CIF] [--tmax TMAX] [--dt DT] [--temp TEMP] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] data`

**Argumente:**

- `data` — Tabelle: V(A^3) E(eV) w1 w2 ... pro Volumen

**Optionen:**

| Option | Beschreibung |
|---|---|
| `-o, --outdir` | Ausgabeordner (Standard: `.`) |
| `--natoms NATOMS` | Atome pro primitiver Zelle, für die Größen pro Atom (Standard: 1) |
| `--cells CELLS` | primitive Zellen pro Superzelle der Moden (Standard: `1`) |
| `--cubic` | zusätzlich a(T). Ohne --structure ist es nur V_prim^(1/3) |
| `--structure CIF` | Struktur des Materials: damit wird a(T) in den KONVENTIONELLEN Gitterparameter umgerechnet (Faktor 4 bei fcc/Diamant, 2 bei bcc) und kubische Symmetrie erkannt |
| `--tmax TMAX` | maximale Temperatur des T-Gitters in K (Standard: 1000) |
| `--dt DT` | Schritt: Integrationsschritt in fs (amorphous) oder Temperaturschritt in K (qha) (Standard: `5.0`) |
| `--temp TEMP` | Arbeitstemperatur in K (Standard: `300.0`) |
| `--dpi DPI` | Auflösung der Bitmap-Abbildungen in Punkten pro Zoll (Standard: 600) |
| `--format` | kommagetrennte Abbildungsformate: pdf,png,svg,eps,tif (Standard: pdf,png) |
| `--no-plot` | nur die Daten exportieren, ohne die Abbildung zu erzeugen |
| `-t, --template` | visuelle Vorlage der Abbildung (die Liste zeigt 'olla-dft templates list') |
| `--font {sans,serif,latex}` | Schriftfamilie: sans, serif oder latex (Computer Modern) |
| `--usetex` | die Texte der Abbildung mit echtem LaTeX setzen (braucht eine LaTeX-Installation) |
| `--palette` | Farbpalette: grayscale, okabe-ito, okabe-ito-dark oder kommagetrennte Hexfarben |
| `--background` | Hintergrundfarbe der Abbildung, z. B. '#FFFFFF' oder 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | Abbildungsbreite je nach Zeitschrift (Standard: generic) |
| `--width` | Breite der Abbildung: single, onehalf oder double, oder eine Zahl in Millimetern |
| `--mono` | Graustufenversion: schwarze Tinte und Linienmuster |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria qha`](THEORY.md)

### `thermochem`

ZPE, Entropie und freie Energie: von einer DFT-Energie zu einer mit dem Experiment vergleichbaren

**Verwendung:** `olla-dft thermochem [-h] [--phase {solido,adsorbato,gas,transicion}] [--structure STRUCTURE] [--temp TEMP] [--pressure PRESSURE] [--symmetry SYMMETRY] [--multiplicity MULTIPLICITY] [--floor FLOOR] [--energy ENERGY] [-o OUTDIR] freqs`

**Argumente:**

- `freqs` — Datei mit Frequenzen in cm-1 oder die kommagetrennte Liste

**Optionen:**

| Option | Beschreibung |
|---|---|
| `--phase {solido,adsorbato,gas,transicion}` | gas fügt Translationen und Rotationen hinzu; transicion verlangt genau eine imaginäre Frequenz (Standard: `solido`) |
| `--structure` | Struktur (nötig für die Gasphase) |
| `--temp TEMP` | Arbeitstemperatur in K (Standard: `298.15`) |
| `--pressure PRESSURE` | in bar (Standard: `1.0`) |
| `--symmetry SYMMETRY` | Symmetriezahl der Punktgruppe: 2 für H2O und O2, 3 für NH3, 12 für CH4 (Standard: `1`) |
| `--multiplicity MULTIPLICITY` | Spinmultiplizität des Grundzustands (Standard: `1`) |
| `--floor FLOOR` | hebt die Moden unterhalb dieses Werts an (cm-1); 100 ist üblich |
| `--energy ENERGY` | E_DFT in eV, um G(T) anzugeben |
| `-o, --outdir` | Ausgabeordner |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria thermochem`](THEORY.md)

### `md`

eine Molekulardynamik-Trajektorie analysieren: g(r), Diffusion und Schwingungsspektrum

**Verwendung:** `olla-dft md [-h] [-o OUTDIR] [--skip SKIP] [--rmax RMAX] [--bins BINS] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] path`

**Argumente:**

- `path` — pw.x-Ausgabe mit calculation='md' oder ihr Ordner

**Optionen:**

| Option | Beschreibung |
|---|---|
| `-o, --outdir` | Ausgabeordner (Standard: `.`) |
| `--skip SKIP` | zu verwerfende Anfangsschritte (Äquilibrierung) |
| `--rmax RMAX` | Cutoff von g(r) in Å; standardmäßig eine halbe Zellkante, bis dahin gilt die Normierung |
| `--bins BINS` | Anzahl der Bins des g(r)-Histogramms (Standard: 200) |
| `--dpi DPI` | Auflösung der Bitmap-Abbildungen in Punkten pro Zoll (Standard: 600) |
| `--format` | kommagetrennte Abbildungsformate: pdf,png,svg,eps,tif (Standard: pdf,png) |
| `--no-plot` | nur die Daten exportieren, ohne die Abbildung zu erzeugen |
| `-t, --template` | visuelle Vorlage der Abbildung (die Liste zeigt 'olla-dft templates list') |
| `--font {sans,serif,latex}` | Schriftfamilie: sans, serif oder latex (Computer Modern) |
| `--usetex` | die Texte der Abbildung mit echtem LaTeX setzen (braucht eine LaTeX-Installation) |
| `--palette` | Farbpalette: grayscale, okabe-ito, okabe-ito-dark oder kommagetrennte Hexfarben |
| `--background` | Hintergrundfarbe der Abbildung, z. B. '#FFFFFF' oder 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | Abbildungsbreite je nach Zeitschrift (Standard: generic) |
| `--width` | Breite der Abbildung: single, onehalf oder double, oder eine Zahl in Millimetern |
| `--mono` | Graustufenversion: schwarze Tinte und Linienmuster |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria md`](THEORY.md)

### `derived`

Debye, Schallgeschwindigkeiten und Slack aus den Cij

**Verwendung:** `olla-dft derived [-h] [--cij CIJ] [--temp TEMP] [-o OUTDIR] file`

**Argumente:**

- `file` — Struktur

**Optionen:**

| Option | Beschreibung |
|---|---|
| `--cij` | Datei mit der elastischen Matrix (Standard: `ELASTIC_C.dat`) |
| `--temp TEMP` | Arbeitstemperatur in K (Standard: `300.0`) |
| `-o, --outdir` | Ordner, in den DERIVED.dat geschrieben wird (Standard: `.`) |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria derived`](THEORY.md)

## Mechanik und Stabilität

### `converge`

Konvergenztests für Cutoffs und k-Gitter

**Verwendung:** `olla-dft converge [-h] [-o OUTDIR] [--run] [--collect] [--pw-cmd PW_CMD] [--nproc NPROC] [-j N] [--redo] [--max-time T] [--estimate] [--timeout TIMEOUT] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--insulator] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--size {paper,poster,presentation}] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--aspect ASPECT] [--mono] [-k {ecutwfc,ecutrho,kmesh}] [--values VALUES] [--threshold THRESHOLD] file`

**Argumente:**

- `file` — Struktur (CIF, POSCAR, pw.x-Eingabe...)

**Optionen:**

| Option | Beschreibung |
|---|---|
| `-o, --outdir` | Ordner des Sweeps (Standard: `convergencia`) |
| `--run` | die Rechnungen jetzt ausführen, eine nach der anderen |
| `--collect` | nur bereits gelaufene Rechnungen analysieren |
| `--pw-cmd` | pw.x-Programm (überschreibt die Konfiguration) |
| `--nproc NPROC` | MPI-Prozesse pro Rechnung |
| `-j, --jobs N` | gleichzeitige Rechnungen (Standard: 1). Ohne --nproc werden die Threads der Maschine unter ihnen aufgeteilt |
| `--redo` | auch bereits abgeschlossene Rechnungen wiederholen |
| `--max-time T` | GESAMTES Zeitbudget: 90m, 2h, 3600. Ist es aufgebraucht, werden keine weiteren gestartet und der Sweep bleibt fortsetzbar |
| `--estimate` | abschätzen, wie lange der Sweep dauern wird, und beenden, anhand der Historie von 'olla-dft db' |
| `--timeout TIMEOUT` | Limit in Sekunden pro Rechnung |
| `--pseudo-dir` | Pseudopotential-Ordner |
| `--pseudo EL=UPF` | ein bestimmtes Pseudopotential erzwingen, z. B. Fe=Fe.rel-pbe.UPF. Kann wiederholt werden. Ohne diese Angabe wählt Olla-DFT mit 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | Wellen-Cutoff (Ry) |
| `--ecutrho ECUTRHO` | Cutoff der Dichte (Ry) |
| `--kspacing KSPACING` | k-Abstand in Å^-1 |
| `--insulator` | occupations='fixed' |
| `--dpi DPI` | Auflösung der Bitmap-Abbildungen in Punkten pro Zoll (Standard: 600) |
| `--format` | kommagetrennte Abbildungsformate: pdf,png,svg,eps,tif (Standard: pdf,png) |
| `--no-plot` | nur die Daten exportieren, ohne die Abbildung zu erzeugen |
| `-t, --template` | visuelle Vorlage der Abbildung |
| `--size {paper,poster,presentation}` | Abbildungsgröße: paper, presentation oder poster |
| `--font {sans,serif,latex}` | Schriftfamilie: sans, serif oder latex (Computer Modern) |
| `--usetex` | die Texte der Abbildung mit echtem LaTeX setzen (braucht eine LaTeX-Installation) |
| `--palette` | Farbpalette: grayscale, okabe-ito, okabe-ito-dark oder kommagetrennte Hexfarben |
| `--background` | Hintergrundfarbe der Abbildung, z. B. '#FFFFFF' oder 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | Abbildungsbreite je nach Zeitschrift (Standard: generic) |
| `--width` | Breite der Abbildung: single, onehalf oder double, oder eine Zahl in Millimetern |
| `--aspect ASPECT` | Verhältnis Höhe/Breite der Abbildung |
| `--mono` | Graustufenversion: schwarze Tinte und Linienmuster |
| `-k, --kind {ecutwfc,ecutrho,kmesh}` | welcher Parameter durchlaufen wird (Standard: ecutwfc) |
| `--values` | kommagetrennte Werte; für kmesh ist 8x8x8 erlaubt |
| `--threshold THRESHOLD` | Konvergenzschwelle in meV/Atom (Standard: 1) |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria converge`](THEORY.md)

### `eos`

E–V-Zustandsgleichung und Kompressionsmodul

**Verwendung:** `olla-dft eos [-h] [-o OUTDIR] [--run] [--collect] [--pw-cmd PW_CMD] [--nproc NPROC] [-j N] [--redo] [--max-time T] [--estimate] [--timeout TIMEOUT] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--insulator] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--size {paper,poster,presentation}] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--aspect ASPECT] [--mono] [--npoints NPOINTS] [--scale SCALE] [--span SPAN] [--equation {birch-murnaghan,murnaghan,vinet}] [--relax-ions] file`

**Argumente:**

- `file` — Struktur (CIF, POSCAR, pw.x-Eingabe...)

**Optionen:**

| Option | Beschreibung |
|---|---|
| `-o, --outdir` | Ordner des Sweeps (Standard: `eos`) |
| `--run` | die Rechnungen jetzt ausführen, eine nach der anderen |
| `--collect` | nur bereits gelaufene Rechnungen analysieren |
| `--pw-cmd` | pw.x-Programm (überschreibt die Konfiguration) |
| `--nproc NPROC` | MPI-Prozesse pro Rechnung |
| `-j, --jobs N` | gleichzeitige Rechnungen (Standard: 1). Ohne --nproc werden die Threads der Maschine unter ihnen aufgeteilt |
| `--redo` | auch bereits abgeschlossene Rechnungen wiederholen |
| `--max-time T` | GESAMTES Zeitbudget: 90m, 2h, 3600. Ist es aufgebraucht, werden keine weiteren gestartet und der Sweep bleibt fortsetzbar |
| `--estimate` | abschätzen, wie lange der Sweep dauern wird, und beenden, anhand der Historie von 'olla-dft db' |
| `--timeout TIMEOUT` | Limit in Sekunden pro Rechnung |
| `--pseudo-dir` | Pseudopotential-Ordner |
| `--pseudo EL=UPF` | ein bestimmtes Pseudopotential erzwingen, z. B. Fe=Fe.rel-pbe.UPF. Kann wiederholt werden. Ohne diese Angabe wählt Olla-DFT mit 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | Wellen-Cutoff (Ry) |
| `--ecutrho ECUTRHO` | Cutoff der Dichte (Ry) |
| `--kspacing KSPACING` | k-Abstand in Å^-1 |
| `--insulator` | occupations='fixed' |
| `--dpi DPI` | Auflösung der Bitmap-Abbildungen in Punkten pro Zoll (Standard: 600) |
| `--format` | kommagetrennte Abbildungsformate: pdf,png,svg,eps,tif (Standard: pdf,png) |
| `--no-plot` | nur die Daten exportieren, ohne die Abbildung zu erzeugen |
| `-t, --template` | visuelle Vorlage der Abbildung |
| `--size {paper,poster,presentation}` | Abbildungsgröße: paper, presentation oder poster |
| `--font {sans,serif,latex}` | Schriftfamilie: sans, serif oder latex (Computer Modern) |
| `--usetex` | die Texte der Abbildung mit echtem LaTeX setzen (braucht eine LaTeX-Installation) |
| `--palette` | Farbpalette: grayscale, okabe-ito, okabe-ito-dark oder kommagetrennte Hexfarben |
| `--background` | Hintergrundfarbe der Abbildung, z. B. '#FFFFFF' oder 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | Abbildungsbreite je nach Zeitschrift (Standard: generic) |
| `--width` | Breite der Abbildung: single, onehalf oder double, oder eine Zahl in Millimetern |
| `--aspect ASPECT` | Verhältnis Höhe/Breite der Abbildung |
| `--mono` | Graustufenversion: schwarze Tinte und Linienmuster |
| `--npoints NPOINTS` | Anzahl der Volumina (Standard: 9) |
| `--scale SCALE` | linearer Faktor, auf den der Sweep zentriert wird (liefert 'olla-dft mlip scan') (Standard: `1.0`) |
| `--span SPAN` | relative Volumenänderung zu jeder Seite (Standard: 0.10) |
| `--equation {birch-murnaghan,murnaghan,vinet}` | Gleichung, die geplottet wird (Standard: `birch-murnaghan`) |
| `--relax-ions` | interne Positionen bei jedem Volumen relaxieren |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria eos`](THEORY.md)

### `elastic`

elastische Konstanten und mechanische Eigenschaften

**Verwendung:** `olla-dft elastic [-h] [-o OUTDIR] [--run] [--collect] [--pw-cmd PW_CMD] [--nproc NPROC] [-j N] [--redo] [--max-time T] [--estimate] [--timeout TIMEOUT] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--insulator] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--size {paper,poster,presentation}] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--aspect ASPECT] [--mono] [--delta DELTA] [--npoints NPOINTS] [--2d] [--thickness A] [--ion-mode {auto,relax,fixed}] file`

**Argumente:**

- `file` — Struktur (CIF, POSCAR, pw.x-Eingabe...)

**Optionen:**

| Option | Beschreibung |
|---|---|
| `-o, --outdir` | Ordner des Sweeps (Standard: `elastic`) |
| `--run` | die Rechnungen jetzt ausführen, eine nach der anderen |
| `--collect` | nur bereits gelaufene Rechnungen analysieren |
| `--pw-cmd` | pw.x-Programm (überschreibt die Konfiguration) |
| `--nproc NPROC` | MPI-Prozesse pro Rechnung |
| `-j, --jobs N` | gleichzeitige Rechnungen (Standard: 1). Ohne --nproc werden die Threads der Maschine unter ihnen aufgeteilt |
| `--redo` | auch bereits abgeschlossene Rechnungen wiederholen |
| `--max-time T` | GESAMTES Zeitbudget: 90m, 2h, 3600. Ist es aufgebraucht, werden keine weiteren gestartet und der Sweep bleibt fortsetzbar |
| `--estimate` | abschätzen, wie lange der Sweep dauern wird, und beenden, anhand der Historie von 'olla-dft db' |
| `--timeout TIMEOUT` | Limit in Sekunden pro Rechnung |
| `--pseudo-dir` | Pseudopotential-Ordner |
| `--pseudo EL=UPF` | ein bestimmtes Pseudopotential erzwingen, z. B. Fe=Fe.rel-pbe.UPF. Kann wiederholt werden. Ohne diese Angabe wählt Olla-DFT mit 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | Wellen-Cutoff (Ry) |
| `--ecutrho ECUTRHO` | Cutoff der Dichte (Ry) |
| `--kspacing KSPACING` | k-Abstand in Å^-1 |
| `--insulator` | occupations='fixed' |
| `--dpi DPI` | Auflösung der Bitmap-Abbildungen in Punkten pro Zoll (Standard: 600) |
| `--format` | kommagetrennte Abbildungsformate: pdf,png,svg,eps,tif (Standard: pdf,png) |
| `--no-plot` | nur die Daten exportieren, ohne die Abbildung zu erzeugen |
| `-t, --template` | visuelle Vorlage der Abbildung |
| `--size {paper,poster,presentation}` | Abbildungsgröße: paper, presentation oder poster |
| `--font {sans,serif,latex}` | Schriftfamilie: sans, serif oder latex (Computer Modern) |
| `--usetex` | die Texte der Abbildung mit echtem LaTeX setzen (braucht eine LaTeX-Installation) |
| `--palette` | Farbpalette: grayscale, okabe-ito, okabe-ito-dark oder kommagetrennte Hexfarben |
| `--background` | Hintergrundfarbe der Abbildung, z. B. '#FFFFFF' oder 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | Abbildungsbreite je nach Zeitschrift (Standard: generic) |
| `--width` | Breite der Abbildung: single, onehalf oder double, oder eine Zahl in Millimetern |
| `--aspect ASPECT` | Verhältnis Höhe/Breite der Abbildung |
| `--mono` | Graustufenversion: schwarze Tinte und Linienmuster |
| `--delta DELTA` | maximal angelegte Dehnung (Standard: 0.010 = 1 %%) |
| `--npoints NPOINTS` | von null verschiedene Dehnungen pro Komponente, gerade Zahl (Standard: 4) |
| `--2d` | Schicht: Konstanten in N/m (nicht in GPa), nur ε1, ε2 und ε6, und Born-Kriterien in 2D |
| `--thickness A` | angenommene Dicke in Å, um auch das GPa-Äquivalent anzugeben (Konvention, keine Messung) |
| `--ion-mode {auto,relax,fixed}` | interne Positionen: auto = fest bei Normaldehnungen und relaxiert bei Scherungen (empfohlen); relax = alle relaxieren; fixed = clamped-ion (Standard: `auto`) |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria elastic`](THEORY.md)

### `strain`

Dehnungs-Sweep: Bandlücke, Energie und Moment als Funktion der angelegten Dehnung

**Verwendung:** `olla-dft strain [-h] [-o OUTDIR] [--run] [--collect] [--pw-cmd PW_CMD] [--nproc NPROC] [-j N] [--redo] [--max-time T] [--estimate] [--timeout TIMEOUT] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--insulator] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--size {paper,poster,presentation}] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--aspect ASPECT] [--mono] [-m {biaxial,cizalla,hidrostatica,uniaxial-a,uniaxial-b,uniaxial-c}] [-r MIN:MAX:N] [--fixed-ions] [--relax-perp] [--nspin {1,2}] [--mag MAG] [--hubbard EL=U] [--vdw {grimme-d2,grimme-d3,DFT-D,ts-vdw,xdm,mbd}] file`

**Argumente:**

- `file` — Struktur (CIF, POSCAR, pw.x-Eingabe...)

**Optionen:**

| Option | Beschreibung |
|---|---|
| `-o, --outdir` | Ordner des Sweeps (Standard: `strain`) |
| `--run` | die Rechnungen jetzt ausführen, eine nach der anderen |
| `--collect` | nur bereits gelaufene Rechnungen analysieren |
| `--pw-cmd` | pw.x-Programm (überschreibt die Konfiguration) |
| `--nproc NPROC` | MPI-Prozesse pro Rechnung |
| `-j, --jobs N` | gleichzeitige Rechnungen (Standard: 1). Ohne --nproc werden die Threads der Maschine unter ihnen aufgeteilt |
| `--redo` | auch bereits abgeschlossene Rechnungen wiederholen |
| `--max-time T` | GESAMTES Zeitbudget: 90m, 2h, 3600. Ist es aufgebraucht, werden keine weiteren gestartet und der Sweep bleibt fortsetzbar |
| `--estimate` | abschätzen, wie lange der Sweep dauern wird, und beenden, anhand der Historie von 'olla-dft db' |
| `--timeout TIMEOUT` | Limit in Sekunden pro Rechnung |
| `--pseudo-dir` | Pseudopotential-Ordner |
| `--pseudo EL=UPF` | ein bestimmtes Pseudopotential erzwingen, z. B. Fe=Fe.rel-pbe.UPF. Kann wiederholt werden. Ohne diese Angabe wählt Olla-DFT mit 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | Wellen-Cutoff (Ry) |
| `--ecutrho ECUTRHO` | Cutoff der Dichte (Ry) |
| `--kspacing KSPACING` | k-Abstand in Å^-1 |
| `--insulator` | occupations='fixed' |
| `--dpi DPI` | Auflösung der Bitmap-Abbildungen in Punkten pro Zoll (Standard: 600) |
| `--format` | kommagetrennte Abbildungsformate: pdf,png,svg,eps,tif (Standard: pdf,png) |
| `--no-plot` | nur die Daten exportieren, ohne die Abbildung zu erzeugen |
| `-t, --template` | visuelle Vorlage der Abbildung |
| `--size {paper,poster,presentation}` | Abbildungsgröße: paper, presentation oder poster |
| `--font {sans,serif,latex}` | Schriftfamilie: sans, serif oder latex (Computer Modern) |
| `--usetex` | die Texte der Abbildung mit echtem LaTeX setzen (braucht eine LaTeX-Installation) |
| `--palette` | Farbpalette: grayscale, okabe-ito, okabe-ito-dark oder kommagetrennte Hexfarben |
| `--background` | Hintergrundfarbe der Abbildung, z. B. '#FFFFFF' oder 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | Abbildungsbreite je nach Zeitschrift (Standard: generic) |
| `--width` | Breite der Abbildung: single, onehalf oder double, oder eine Zahl in Millimetern |
| `--aspect ASPECT` | Verhältnis Höhe/Breite der Abbildung |
| `--mono` | Graustufenversion: schwarze Tinte und Linienmuster |
| `-m, --mode {biaxial,cizalla,hidrostatica,uniaxial-a,uniaxial-b,uniaxial-c}` | was gedehnt wird (Standard: biaxial) |
| `-r, --range MIN:MAX:N` | Bereich in PROZENT, z. B. -5:5:11 (von -5 %% bis +5 %% in 11 Punkten) (Standard: `-5:5:11`) |
| `--fixed-ions` | die internen Positionen bei jeder Dehnung nicht relaxieren (schneller und weniger realistisch) |
| `--relax-perp` | die Achse senkrecht zur gedehnten Ebene freilassen (Poisson-Relaxation); bei Schichten unverzichtbar |
| `--nspin {1,2}` | 2 schaltet die Spinpolarisation ein (Standard: `1`) |
| `--mag` | Startmagnetisierung (impliziert --nspin 2) |
| `--hubbard EL=U` | Hubbard-U in eV pro Element |
| `--vdw {grimme-d2,grimme-d3,DFT-D,ts-vdw,xdm,mbd}` | Dispersionskorrektur |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria strain`](THEORY.md)

### `layers`

Lagen, basalen Abstand und Zwischenschichtlücke erkennen

**Verwendung:** `olla-dft layers [-h] [--tol TOL] [--wavelength WAVELENGTH] [--slab ARCHIVO] [--vacuum VACUUM] file`

**Argumente:**

- `file` — Eingabestruktur (CIF, POSCAR, pw.x-Eingabe...)

**Optionen:**

| Option | Beschreibung |
|---|---|
| `--tol TOL` | Bindungstoleranz zusätzlich zu den Kovalenzradien (Å) (Standard: `0.45`) |
| `--wavelength` | Strahlung für die basalen Reflexe (Standard CuKa) |
| `--slab ARCHIVO` | zusätzlich die Monolage mit Vakuum in diese Datei schreiben |
| `--vacuum VACUUM` | Vakuum der Monolage in Å (Standard 20) |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria layers`](THEORY.md)

### `xrd`

simuliertes Pulverdiffraktogramm

**Verwendung:** `olla-dft xrd [-h] [-o OUTDIR] [--suite] [--basis {conventional,input}] [--wavelength WAVELENGTH] [--tt-min TT_MIN] [--tt-max TT_MAX] [--fwhm FWHM] [--size SIZE] [--biso BISO] [--exp EXP] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--size-preset {paper,poster,presentation}] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--aspect ASPECT] [--mono] file`

**Argumente:**

- `file` — Eingabestruktur (CIF, POSCAR, pw.x-Eingabe...)

**Optionen:**

| Option | Beschreibung |
|---|---|
| `-o, --outdir` | Ausgabeordner (Standard: `.`) |
| `--suite` | zusätzlich ein Austausch-JSON für die anderen Apps der Suite exportieren |
| `--basis {conventional,input}` | Zelle, in der die hkl indiziert werden: 'conventional' (Standard, die Indizes der PDF-Karten) oder 'input' (die Zelle der Datei, wie sie ist) |
| `--wavelength` | Strahlung: AgKa, CoKa, CrKa, CuKa, CuKa1, FeKa, MoKa oder λ in Å (Standard: `CuKa`) |
| `--tt-min TT_MIN` | minimales 2θ (°) (Standard: `5.0`) |
| `--tt-max TT_MAX` | maximales 2θ (°) (Standard: `70.0`) |
| `--fwhm FWHM` | instrumentelle Breite (° 2θ, Standard 0.15) |
| `--size SIZE` | Kristallitgröße in nm (Scherrer-Verbreiterung) |
| `--biso BISO` | globaler Temperaturfaktor B (Å²) |
| `--exp` | experimentelles Diffraktogramm (2θ, I) zum Überlagern |
| `--dpi DPI` | Auflösung der Bitmap-Abbildungen in Punkten pro Zoll (Standard: 600) |
| `--format` | kommagetrennte Abbildungsformate: pdf,png,svg,eps,tif (Standard: pdf,png) |
| `--no-plot` | nur die Daten exportieren, ohne die Abbildung zu erzeugen |
| `-t, --template` | visuelle Vorlage der Abbildung (die Liste zeigt 'olla-dft templates list') |
| `--size-preset {paper,poster,presentation}` | Schriftskala der Abbildung |
| `--font {sans,serif,latex}` | Schriftfamilie: sans, serif oder latex (Computer Modern) |
| `--usetex` | die Texte der Abbildung mit echtem LaTeX setzen (braucht eine LaTeX-Installation) |
| `--palette` | Farbpalette: grayscale, okabe-ito, okabe-ito-dark oder kommagetrennte Hexfarben |
| `--background` | Hintergrundfarbe der Abbildung, z. B. '#FFFFFF' oder 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | Abbildungsbreite je nach Zeitschrift (Standard: generic) |
| `--width` | Breite der Abbildung: single, onehalf oder double, oder eine Zahl in Millimetern |
| `--aspect ASPECT` | Verhältnis Höhe/Breite der Abbildung |
| `--mono` | Graustufenversion: schwarze Tinte und Linienmuster |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria xrd`](THEORY.md)

### `exfoliate`

Exfoliationsenergie (Bulk vs Monolage)

**Verwendung:** `olla-dft exfoliate [-h] [-o OUTDIR] [--run] [--collect] [--pw-cmd PW_CMD] [--nproc NPROC] [-j N] [--redo] [--max-time T] [--estimate] [--timeout TIMEOUT] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--insulator] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--size {paper,poster,presentation}] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--aspect ASPECT] [--mono] [--vacuum VACUUM] [--vdw {grimme-d2,grimme-d3,DFT-D,ts-vdw,xdm,mbd}] [--tol TOL] [--relax-slab] file`

**Argumente:**

- `file` — Struktur (CIF, POSCAR, pw.x-Eingabe...)

**Optionen:**

| Option | Beschreibung |
|---|---|
| `-o, --outdir` | Ordner des Sweeps (Standard: `exfoliacion`) |
| `--run` | die Rechnungen jetzt ausführen, eine nach der anderen |
| `--collect` | nur bereits gelaufene Rechnungen analysieren |
| `--pw-cmd` | pw.x-Programm (überschreibt die Konfiguration) |
| `--nproc NPROC` | MPI-Prozesse pro Rechnung |
| `-j, --jobs N` | gleichzeitige Rechnungen (Standard: 1). Ohne --nproc werden die Threads der Maschine unter ihnen aufgeteilt |
| `--redo` | auch bereits abgeschlossene Rechnungen wiederholen |
| `--max-time T` | GESAMTES Zeitbudget: 90m, 2h, 3600. Ist es aufgebraucht, werden keine weiteren gestartet und der Sweep bleibt fortsetzbar |
| `--estimate` | abschätzen, wie lange der Sweep dauern wird, und beenden, anhand der Historie von 'olla-dft db' |
| `--timeout TIMEOUT` | Limit in Sekunden pro Rechnung |
| `--pseudo-dir` | Pseudopotential-Ordner |
| `--pseudo EL=UPF` | ein bestimmtes Pseudopotential erzwingen, z. B. Fe=Fe.rel-pbe.UPF. Kann wiederholt werden. Ohne diese Angabe wählt Olla-DFT mit 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | Wellen-Cutoff (Ry) |
| `--ecutrho ECUTRHO` | Cutoff der Dichte (Ry) |
| `--kspacing KSPACING` | k-Abstand in Å^-1 |
| `--insulator` | occupations='fixed' |
| `--dpi DPI` | Auflösung der Bitmap-Abbildungen in Punkten pro Zoll (Standard: 600) |
| `--format` | kommagetrennte Abbildungsformate: pdf,png,svg,eps,tif (Standard: pdf,png) |
| `--no-plot` | nur die Daten exportieren, ohne die Abbildung zu erzeugen |
| `-t, --template` | visuelle Vorlage der Abbildung |
| `--size {paper,poster,presentation}` | Abbildungsgröße: paper, presentation oder poster |
| `--font {sans,serif,latex}` | Schriftfamilie: sans, serif oder latex (Computer Modern) |
| `--usetex` | die Texte der Abbildung mit echtem LaTeX setzen (braucht eine LaTeX-Installation) |
| `--palette` | Farbpalette: grayscale, okabe-ito, okabe-ito-dark oder kommagetrennte Hexfarben |
| `--background` | Hintergrundfarbe der Abbildung, z. B. '#FFFFFF' oder 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | Abbildungsbreite je nach Zeitschrift (Standard: generic) |
| `--width` | Breite der Abbildung: single, onehalf oder double, oder eine Zahl in Millimetern |
| `--aspect ASPECT` | Verhältnis Höhe/Breite der Abbildung |
| `--mono` | Graustufenversion: schwarze Tinte und Linienmuster |
| `--vacuum VACUUM` | Vakuum der Monolage in Å (Standard 20) |
| `--vdw {grimme-d2,grimme-d3,DFT-D,ts-vdw,xdm,mbd}` | Dispersionskorrektur für beide Rechnungen |
| `--tol TOL` | Bindungstoleranz zur Erkennung der Lagen (Å) (Standard: `0.45`) |
| `--relax-slab` | die Positionen der Monolage relaxieren |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria exfoliate`](THEORY.md)

### `gamma`

Oberflächen- und Spaltenergie über den linearen Fit nach Fiorentini–Methfessel, mit Konvergenz gegen die Slab-Dicke

**Verwendung:** `olla-dft gamma [-h] [-o OUTDIR] [--run] [--collect] [--pw-cmd PW_CMD] [--nproc NPROC] [-j N] [--redo] [--max-time T] [--estimate] [--timeout TIMEOUT] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--insulator] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--size {paper,poster,presentation}] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--aspect ASPECT] [--mono] [-m MILLER] [-l LAYERS] [--vacuum VACUUM] [--fix N] [--relax] [--no-bulk] [--no-reduce] [--vdw {grimme-d2,grimme-d3,DFT-D,ts-vdw,xdm,mbd}] [--dipole] [--nspin {1,2}] [--mag MAG] file`

**Argumente:**

- `file` — Struktur (CIF, POSCAR, pw.x-Eingabe...)

**Optionen:**

| Option | Beschreibung |
|---|---|
| `-o, --outdir` | Ordner des Sweeps (Standard: `gamma`) |
| `--run` | die Rechnungen jetzt ausführen, eine nach der anderen |
| `--collect` | nur bereits gelaufene Rechnungen analysieren |
| `--pw-cmd` | pw.x-Programm (überschreibt die Konfiguration) |
| `--nproc NPROC` | MPI-Prozesse pro Rechnung |
| `-j, --jobs N` | gleichzeitige Rechnungen (Standard: 1). Ohne --nproc werden die Threads der Maschine unter ihnen aufgeteilt |
| `--redo` | auch bereits abgeschlossene Rechnungen wiederholen |
| `--max-time T` | GESAMTES Zeitbudget: 90m, 2h, 3600. Ist es aufgebraucht, werden keine weiteren gestartet und der Sweep bleibt fortsetzbar |
| `--estimate` | abschätzen, wie lange der Sweep dauern wird, und beenden, anhand der Historie von 'olla-dft db' |
| `--timeout TIMEOUT` | Limit in Sekunden pro Rechnung |
| `--pseudo-dir` | Pseudopotential-Ordner |
| `--pseudo EL=UPF` | ein bestimmtes Pseudopotential erzwingen, z. B. Fe=Fe.rel-pbe.UPF. Kann wiederholt werden. Ohne diese Angabe wählt Olla-DFT mit 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | Wellen-Cutoff (Ry) |
| `--ecutrho ECUTRHO` | Cutoff der Dichte (Ry) |
| `--kspacing KSPACING` | k-Abstand in Å^-1 |
| `--insulator` | occupations='fixed' |
| `--dpi DPI` | Auflösung der Bitmap-Abbildungen in Punkten pro Zoll (Standard: 600) |
| `--format` | kommagetrennte Abbildungsformate: pdf,png,svg,eps,tif (Standard: pdf,png) |
| `--no-plot` | nur die Daten exportieren, ohne die Abbildung zu erzeugen |
| `-t, --template` | visuelle Vorlage der Abbildung |
| `--size {paper,poster,presentation}` | Abbildungsgröße: paper, presentation oder poster |
| `--font {sans,serif,latex}` | Schriftfamilie: sans, serif oder latex (Computer Modern) |
| `--usetex` | die Texte der Abbildung mit echtem LaTeX setzen (braucht eine LaTeX-Installation) |
| `--palette` | Farbpalette: grayscale, okabe-ito, okabe-ito-dark oder kommagetrennte Hexfarben |
| `--background` | Hintergrundfarbe der Abbildung, z. B. '#FFFFFF' oder 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | Abbildungsbreite je nach Zeitschrift (Standard: generic) |
| `--width` | Breite der Abbildung: single, onehalf oder double, oder eine Zahl in Millimetern |
| `--aspect ASPECT` | Verhältnis Höhe/Breite der Abbildung |
| `--mono` | Graustufenversion: schwarze Tinte und Linienmuster |
| `-m, --miller` | Miller-Indizes der Facette, z. B. '1 1 1' (Standard: `1 0 0`) |
| `-l, --layers` | zu berechnende Dicken, kommagetrennt (Standard: 3,4,5,6). Mindestens zwei sind nötig |
| `--vacuum VACUUM` | Vakuum in Å (Standard: 20) |
| `--fix N` | beim Relaxieren N untere Lagen einfrieren |
| `--relax` | die Positionen relaxieren (γ sinkt um 5 bis 20 %%) |
| `--no-bulk` | das Bulk nicht separat berechnen; nur der lineare Fit E_slab(N) = 2γA + N·E_bulk |
| `--no-reduce` | die Oberflächenzelle nicht auf die minimale reduzieren (standardmäßig wird reduziert: gleiches γ, viel billiger) |
| `--vdw {grimme-d2,grimme-d3,DFT-D,ts-vdw,xdm,mbd}` | Dispersionskorrektur (van der Waals): grimme-d2, grimme-d3, DFT-D, ts-vdw, xdm oder mbd |
| `--dipole` | Dipolkorrektur, für polare Slabs |
| `--nspin {1,2}` | 2 schaltet die Spinpolarisation ein (Standard: 1) |
| `--mag` | Startmagnetisierung (impliziert --nspin 2) |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria gamma`](THEORY.md)

## Oberflächen, Defekte und Chemie

### `surface`

eine (hkl)-Oberfläche mit Vakuum schneiden

**Verwendung:** `olla-dft surface [-h] [-m MILLER] [-l LAYERS] [--vacuum VACUUM] [--fix FIX] [-o OUTPUT] file`

**Argumente:**

- `file` — Eingabestruktur (CIF, POSCAR, pw.x-Eingabe...)

**Optionen:**

| Option | Beschreibung |
|---|---|
| `-m, --miller` | Miller-Indizes, z. B. '1 1 1' oder 1,1,1 (Standard: `1 0 0`) |
| `-l, --layers LAYERS` | Anzahl der Atomlagen des Slabs (Standard: 6) |
| `--vacuum VACUUM` | Gesamtvakuum in Å (Standard: `15.0`) |
| `--fix FIX` | einzufrierende untere Atomebenen |
| `-o, --output` | Ausgabedatei (CIF/POSCAR) |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria surface`](THEORY.md)

### `defect`

einen Punktdefekt erzeugen

**Verwendung:** `olla-dft defect [-h] [-k {vacancy,substitution,interstitial}] [--site SITE] [--new-element NEW_ELEMENT] [--supercell SUPERCELL] [--position POSITION] [-o OUTDIR] file`

**Argumente:**

- `file` — Eingabestruktur (CIF, POSCAR, pw.x-Eingabe...)

**Optionen:**

| Option | Beschreibung |
|---|---|
| `-k, --kind {vacancy,substitution,interstitial}` | Defekttyp: vacancy, substitution oder interstitial (Standard: vacancy) |
| `--site SITE` | Index des betroffenen Atoms (ab 0) |
| `--new-element` | eintretende Spezies |
| `--supercell` | z. B. 3x3x3 |
| `--position` | fraktionelle x,y,z (interstitiell) |
| `-o, --outdir` | Ausgabeordner (Standard: `defecto`) |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria defect`](THEORY.md)

### `interface`

Heterostruktur: zwei Materialien mit der kleinstmöglichen Gitterverzerrung stapeln

**Verwendung:** `olla-dft interface [-h] [-o OUTDIR] [--name NAME] [--max-index MAX_INDEX] [--tol TOL] [--max-atoms MAX_ATOMS] [--index INDEX] [--top TOP] [--list] [--separation SEPARATION] [--vacuum VACUUM] [--strain {first,second,both}] [--shift SHIFT] file1 file2`

**Argumente:**

- `file1` — unteres Material (das Substrat)
- `file2` — oberes Material

**Optionen:**

| Option | Beschreibung |
|---|---|
| `-o, --outdir` | Ausgabeordner (Standard: `.`) |
| `--name` | Basisname der Heterostruktur-Dateien (Standard: heteroestructura) |
| `--max-index MAX_INDEX` | größter ganzzahliger Koeffizient der Superzelle; ihn zu erhöhen findet stärker verdrehte Zellen, dauert aber viel länger (Standard: `4`) |
| `--tol TOL` | maximal akzeptierte Dehnung (0.05 = 5 %%) (Standard: `0.05`) |
| `--max-atoms MAX_ATOMS` | maximal zulässige Atomzahl in der Superzelle der Grenzfläche (Standard: 200) |
| `--index INDEX` | welche der Kandidaten gebaut wird |
| `--top TOP` | wie viele Kandidaten aufgelistet werden, von kleinster zu größter Dehnung (Standard: 10) |
| `--list` | nur die Kandidaten auflisten, ohne etwas zu bauen |
| `--separation SEPARATION` | Anfangsabstand zwischen den Lagen in Å; standardmäßig aus den Van-der-Waals-Radien |
| `--vacuum VACUUM` | Vakuum über der Heterostruktur in Å (Standard: 20) |
| `--strain {first,second,both}` | wer gedehnt wird: das untere, das obere oder beide je zur Hälfte (Standard: `second`) |
| `--shift` | seitliche Verschiebung des oberen Materials, in Bruchteilen der gemeinsamen Zelle |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria interface`](THEORY.md)

### `adsorb`

Adsorptionsplätze auf einem Slab und ihre Energie

**Verwendung:** `olla-dft adsorb [-h] [-o OUTDIR] [--run] [--collect] [--pw-cmd PW_CMD] [--nproc NPROC] [-j N] [--redo] [--max-time T] [--estimate] [--timeout TIMEOUT] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--insulator] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--size {paper,poster,presentation}] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--aspect ASPECT] [--mono] --mol MOLECULA [--sites SITES] [--height HEIGHT] [--face {top,bottom}] [--rotations ROTATIONS] [--anchor ANCHOR] [--fixed-ions] [--vdw {grimme-d2,grimme-d3,DFT-D,ts-vdw,xdm,mbd}] [--dipole] [--nspin {1,2}] [--mag MAG] file`

**Argumente:**

- `file` — Struktur (CIF, POSCAR, pw.x-Eingabe...)

**Optionen:**

| Option | Beschreibung |
|---|---|
| `-o, --outdir` | Ordner des Sweeps (Standard: `adsorb`) |
| `--run` | die Rechnungen jetzt ausführen, eine nach der anderen |
| `--collect` | nur bereits gelaufene Rechnungen analysieren |
| `--pw-cmd` | pw.x-Programm (überschreibt die Konfiguration) |
| `--nproc NPROC` | MPI-Prozesse pro Rechnung |
| `-j, --jobs N` | gleichzeitige Rechnungen (Standard: 1). Ohne --nproc werden die Threads der Maschine unter ihnen aufgeteilt |
| `--redo` | auch bereits abgeschlossene Rechnungen wiederholen |
| `--max-time T` | GESAMTES Zeitbudget: 90m, 2h, 3600. Ist es aufgebraucht, werden keine weiteren gestartet und der Sweep bleibt fortsetzbar |
| `--estimate` | abschätzen, wie lange der Sweep dauern wird, und beenden, anhand der Historie von 'olla-dft db' |
| `--timeout TIMEOUT` | Limit in Sekunden pro Rechnung |
| `--pseudo-dir` | Pseudopotential-Ordner |
| `--pseudo EL=UPF` | ein bestimmtes Pseudopotential erzwingen, z. B. Fe=Fe.rel-pbe.UPF. Kann wiederholt werden. Ohne diese Angabe wählt Olla-DFT mit 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | Wellen-Cutoff (Ry) |
| `--ecutrho ECUTRHO` | Cutoff der Dichte (Ry) |
| `--kspacing KSPACING` | k-Abstand in Å^-1 |
| `--insulator` | occupations='fixed' |
| `--dpi DPI` | Auflösung der Bitmap-Abbildungen in Punkten pro Zoll (Standard: 600) |
| `--format` | kommagetrennte Abbildungsformate: pdf,png,svg,eps,tif (Standard: pdf,png) |
| `--no-plot` | nur die Daten exportieren, ohne die Abbildung zu erzeugen |
| `-t, --template` | visuelle Vorlage der Abbildung |
| `--size {paper,poster,presentation}` | Abbildungsgröße: paper, presentation oder poster |
| `--font {sans,serif,latex}` | Schriftfamilie: sans, serif oder latex (Computer Modern) |
| `--usetex` | die Texte der Abbildung mit echtem LaTeX setzen (braucht eine LaTeX-Installation) |
| `--palette` | Farbpalette: grayscale, okabe-ito, okabe-ito-dark oder kommagetrennte Hexfarben |
| `--background` | Hintergrundfarbe der Abbildung, z. B. '#FFFFFF' oder 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | Abbildungsbreite je nach Zeitschrift (Standard: generic) |
| `--width` | Breite der Abbildung: single, onehalf oder double, oder eine Zahl in Millimetern |
| `--aspect ASPECT` | Verhältnis Höhe/Breite der Abbildung |
| `--mono` | Graustufenversion: schwarze Tinte und Linienmuster |
| `--mol MOLECULA` | Adsorbat: ein Name aus der ASE-Datenbank (CO2, H2O, CO, NH3, O2...) oder eine Datei mit dem Molekül |
| `--sites` | zu testende Platztypen (Standard: alle drei) |
| `--height HEIGHT` | Anfangshöhe des Adsorbats über dem Platz, in Å (Standard: 2.0) |
| `--face {top,bottom}` | Seite des Slabs, auf der adsorbiert wird (Standard: `top`) |
| `--rotations ROTATIONS` | zu testende Orientierungen durch Drehung um die Normale (Standard: 1) |
| `--anchor ANCHOR` | Atom des Moleküls, das auf dem Platz sitzt (Index ab 0; Standard: 0) |
| `--fixed-ions` | nicht relaxieren: nur scf an der Anfangsgeometrie |
| `--vdw {grimme-d2,grimme-d3,DFT-D,ts-vdw,xdm,mbd}` | Dispersionskorrektur (bei Physisorption fast Pflicht) |
| `--dipole` | Dipolkorrektur: ein Slab mit Adsorbat auf nur einer Seite ist polar |
| `--nspin {1,2}` | 2 schaltet die Spinpolarisation ein (Standard: 1) |
| `--mag` | Startmagnetisierung (impliziert --nspin 2) |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria adsorb`](THEORY.md)

### `eform`

Bildungsenergie geladener Defekte, Übergangsniveaus und Diagramm E_f vs ε_F

**Verwendung:** `olla-dft eform [-h] [-o OUTDIR] [--run] [--collect] [--pw-cmd PW_CMD] [--nproc NPROC] [-j N] [--redo] [--max-time T] [--estimate] [--timeout TIMEOUT] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--insulator] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--size {paper,poster,presentation}] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--aspect ASPECT] [--mono] [-k {vacancy,substitution,interstitial}] [--site SITE] [--new-element NEW_ELEMENT] [--position POSITION] [--supercell SUPERCELL] [-q CHARGES] [--epsilon EPSILON] [--correction {ninguna,makov-payne,lany-zunger}] [--mu EL=eV] [--align POT_DEF POT_PERF] [--dv DV] [--fixed-ions] [--vdw {grimme-d2,grimme-d3,DFT-D,ts-vdw,xdm,mbd}] [--nspin {1,2}] [--mag MAG] file`

**Argumente:**

- `file` — Struktur (CIF, POSCAR, pw.x-Eingabe...)

**Optionen:**

| Option | Beschreibung |
|---|---|
| `-o, --outdir` | Ordner des Sweeps (Standard: `formacion`) |
| `--run` | die Rechnungen jetzt ausführen, eine nach der anderen |
| `--collect` | nur bereits gelaufene Rechnungen analysieren |
| `--pw-cmd` | pw.x-Programm (überschreibt die Konfiguration) |
| `--nproc NPROC` | MPI-Prozesse pro Rechnung |
| `-j, --jobs N` | gleichzeitige Rechnungen (Standard: 1). Ohne --nproc werden die Threads der Maschine unter ihnen aufgeteilt |
| `--redo` | auch bereits abgeschlossene Rechnungen wiederholen |
| `--max-time T` | GESAMTES Zeitbudget: 90m, 2h, 3600. Ist es aufgebraucht, werden keine weiteren gestartet und der Sweep bleibt fortsetzbar |
| `--estimate` | abschätzen, wie lange der Sweep dauern wird, und beenden, anhand der Historie von 'olla-dft db' |
| `--timeout TIMEOUT` | Limit in Sekunden pro Rechnung |
| `--pseudo-dir` | Pseudopotential-Ordner |
| `--pseudo EL=UPF` | ein bestimmtes Pseudopotential erzwingen, z. B. Fe=Fe.rel-pbe.UPF. Kann wiederholt werden. Ohne diese Angabe wählt Olla-DFT mit 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | Wellen-Cutoff (Ry) |
| `--ecutrho ECUTRHO` | Cutoff der Dichte (Ry) |
| `--kspacing KSPACING` | k-Abstand in Å^-1 |
| `--insulator` | occupations='fixed' |
| `--dpi DPI` | Auflösung der Bitmap-Abbildungen in Punkten pro Zoll (Standard: 600) |
| `--format` | kommagetrennte Abbildungsformate: pdf,png,svg,eps,tif (Standard: pdf,png) |
| `--no-plot` | nur die Daten exportieren, ohne die Abbildung zu erzeugen |
| `-t, --template` | visuelle Vorlage der Abbildung |
| `--size {paper,poster,presentation}` | Abbildungsgröße: paper, presentation oder poster |
| `--font {sans,serif,latex}` | Schriftfamilie: sans, serif oder latex (Computer Modern) |
| `--usetex` | die Texte der Abbildung mit echtem LaTeX setzen (braucht eine LaTeX-Installation) |
| `--palette` | Farbpalette: grayscale, okabe-ito, okabe-ito-dark oder kommagetrennte Hexfarben |
| `--background` | Hintergrundfarbe der Abbildung, z. B. '#FFFFFF' oder 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | Abbildungsbreite je nach Zeitschrift (Standard: generic) |
| `--width` | Breite der Abbildung: single, onehalf oder double, oder eine Zahl in Millimetern |
| `--aspect ASPECT` | Verhältnis Höhe/Breite der Abbildung |
| `--mono` | Graustufenversion: schwarze Tinte und Linienmuster |
| `-k, --kind {vacancy,substitution,interstitial}` | Defekttyp (Standard: `vacancy`) |
| `--site SITE` | Index des betroffenen Atoms in der Superzelle (ab 0) |
| `--new-element` | eintretende Spezies |
| `--position` | fraktionelle x,y,z (interstitiell) |
| `--supercell` | Größe der Superzelle (Standard: 2x2x2) |
| `-q, --charges` | kommagetrennte Ladungszustände, z. B. -2,-1,0,1,2 (Standard: `0`) |
| `--epsilon EPSILON` | Dielektrizitätskonstante des Materials, zur Abschirmung der Bildladungskorrektur |
| `--correction {ninguna,makov-payne,lany-zunger}` | Schema der Finite-Size-Korrektur (Standard: `lany-zunger`) |
| `--mu EL=eV` | chemisches Potential pro Element, in eV pro Atom. Kann wiederholt werden |
| `--align ('POT_DEF', 'POT_PERF')` | zwei Cube-Dateien des elektrostatischen Potentials (defekt und perfekt) für den ΔV-Term |
| `--dv DV` | ΔV-Ausrichtung in eV, falls du sie schon berechnet hast |
| `--fixed-ions` | den Defekt nicht in jedem Ladungszustand relaxieren |
| `--vdw {grimme-d2,grimme-d3,DFT-D,ts-vdw,xdm,mbd}` | Dispersionskorrektur (van der Waals): grimme-d2, grimme-d3, DFT-D, ts-vdw, xdm oder mbd |
| `--nspin {1,2}` | 2 schaltet die Spinpolarisation ein (Standard: 1) |
| `--mag` | Startmagnetisierung (impliziert --nspin 2) |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria eform`](THEORY.md)

### `align`

Bandausrichtung zwischen zwei Materialien: Offsets ΔE_v, ΔE_c und Typ I/II/III

**Verwendung:** `olla-dft align [-h] [--interface CARPETA] [--names NAMES] [--axis AXIS] [--window A] [--rerun] [-o OUTDIR] [--pw-cmd PW_CMD] [--nproc NPROC] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] a b`

**Argumente:**

- `a` — Rechnungsordner des ersten Materials
- `b` — Rechnungsordner des zweiten Materials

**Optionen:**

| Option | Beschreibung |
|---|---|
| `--interface CARPETA` | Ordner der Grenzfläche; aktiviert die rigorose Methode von Van de Walle-Martin |
| `--names` | Namen für den Bericht, kommagetrennt (standardmäßig die Ordnernamen) |
| `--axis` | Achse des planaren Profils (Standard: `c`) |
| `--window A` | Fenster des makroskopischen Mittels in Å (standardmäßig ein Achtel der Zelle) |
| `--rerun` | pp.x erneut ausführen, auch wenn die Cube-Datei bereits existiert |
| `-o, --outdir` | Ausgabeordner (Standard: `alineamiento`) |
| `--pw-cmd` | pw.x-Programm für --run; aus seinem Pfad werden die übrigen QE-Binärdateien gefunden |
| `--nproc NPROC` | Anzahl der MPI-Prozesse für die mit --run gestarteten Rechnungen |
| `--dpi DPI` | Auflösung der Bitmap-Abbildungen in Punkten pro Zoll (Standard: 600) |
| `--format` | kommagetrennte Abbildungsformate: pdf,png,svg,eps,tif (Standard: pdf,png) |
| `--no-plot` | nur die Daten exportieren, ohne die Abbildung zu erzeugen |
| `-t, --template` | visuelle Vorlage der Abbildung (die Liste zeigt 'olla-dft templates list') |
| `--font {sans,serif,latex}` | Schriftfamilie: sans, serif oder latex (Computer Modern) |
| `--usetex` | die Texte der Abbildung mit echtem LaTeX setzen (braucht eine LaTeX-Installation) |
| `--palette` | Farbpalette: grayscale, okabe-ito, okabe-ito-dark oder kommagetrennte Hexfarben |
| `--background` | Hintergrundfarbe der Abbildung, z. B. '#FFFFFF' oder 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | Abbildungsbreite je nach Zeitschrift (Standard: generic) |
| `--width` | Breite der Abbildung: single, onehalf oder double, oder eine Zahl in Millimetern |
| `--mono` | Graustufenversion: schwarze Tinte und Linienmuster |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria align`](THEORY.md)

### `esm`

geladene Oberflächen mit effektivem Abschirmmedium: Austrittsarbeit, Kapazität und Potential der Nullladung

**Verwendung:** `olla-dft esm [-h] [-o OUTDIR] [--run] [--collect] [--pw-cmd PW_CMD] [--nproc NPROC] [-j N] [--redo] [--max-time T] [--estimate] [--timeout TIMEOUT] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--insulator] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--size {paper,poster,presentation}] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--aspect ASPECT] [--mono] [--bc {bc1,bc2,bc3}] [--charge CHARGE] [--field FIELD] [--esm-w WIDTH_ESM] [--nfit NFIT] file`

**Argumente:**

- `file` — Struktur (CIF, POSCAR, pw.x-Eingabe...)

**Optionen:**

| Option | Beschreibung |
|---|---|
| `-o, --outdir` | Ordner des Sweeps (Standard: `esm`) |
| `--run` | die Rechnungen jetzt ausführen, eine nach der anderen |
| `--collect` | nur bereits gelaufene Rechnungen analysieren |
| `--pw-cmd` | pw.x-Programm (überschreibt die Konfiguration) |
| `--nproc NPROC` | MPI-Prozesse pro Rechnung |
| `-j, --jobs N` | gleichzeitige Rechnungen (Standard: 1). Ohne --nproc werden die Threads der Maschine unter ihnen aufgeteilt |
| `--redo` | auch bereits abgeschlossene Rechnungen wiederholen |
| `--max-time T` | GESAMTES Zeitbudget: 90m, 2h, 3600. Ist es aufgebraucht, werden keine weiteren gestartet und der Sweep bleibt fortsetzbar |
| `--estimate` | abschätzen, wie lange der Sweep dauern wird, und beenden, anhand der Historie von 'olla-dft db' |
| `--timeout TIMEOUT` | Limit in Sekunden pro Rechnung |
| `--pseudo-dir` | Pseudopotential-Ordner |
| `--pseudo EL=UPF` | ein bestimmtes Pseudopotential erzwingen, z. B. Fe=Fe.rel-pbe.UPF. Kann wiederholt werden. Ohne diese Angabe wählt Olla-DFT mit 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | Wellen-Cutoff (Ry) |
| `--ecutrho ECUTRHO` | Cutoff der Dichte (Ry) |
| `--kspacing KSPACING` | k-Abstand in Å^-1 |
| `--insulator` | occupations='fixed' |
| `--dpi DPI` | Auflösung der Bitmap-Abbildungen in Punkten pro Zoll (Standard: 600) |
| `--format` | kommagetrennte Abbildungsformate: pdf,png,svg,eps,tif (Standard: pdf,png) |
| `--no-plot` | nur die Daten exportieren, ohne die Abbildung zu erzeugen |
| `-t, --template` | visuelle Vorlage der Abbildung |
| `--size {paper,poster,presentation}` | Abbildungsgröße: paper, presentation oder poster |
| `--font {sans,serif,latex}` | Schriftfamilie: sans, serif oder latex (Computer Modern) |
| `--usetex` | die Texte der Abbildung mit echtem LaTeX setzen (braucht eine LaTeX-Installation) |
| `--palette` | Farbpalette: grayscale, okabe-ito, okabe-ito-dark oder kommagetrennte Hexfarben |
| `--background` | Hintergrundfarbe der Abbildung, z. B. '#FFFFFF' oder 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | Abbildungsbreite je nach Zeitschrift (Standard: generic) |
| `--width` | Breite der Abbildung: single, onehalf oder double, oder eine Zahl in Millimetern |
| `--aspect ASPECT` | Verhältnis Höhe/Breite der Abbildung |
| `--mono` | Graustufenversion: schwarze Tinte und Linienmuster |
| `--bc {bc1,bc2,bc3}` | bc1 Vakuum/Vakuum (neutrale Slabs), bc2 Metall/Metall (Kondensator), bc3 Vakuum/Metall (Elektrode, neben bc2 die einzige, die eine Nettoladung erlaubt) (Standard: `bc1`) |
| `--charge` | Nettoladungen in e, kommagetrennt: -0.2,0,0.2 (Standard: `0`) |
| `--field FIELD` | angelegtes Feld in Ry/a.u. (nur mit bc2) |
| `--esm-w WIDTH_ESM` | Verschiebung der ESM-Grenze in a.u. |
| `--nfit NFIT` | Fitpunkte des Potentials an der Grenze (Standard: `4`) |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria esm`](THEORY.md)

### `echem`

computergestützte Wasserstoffelektrode: HER, OER, limitierendes Potential und Überpotential

**Verwendung:** `olla-dft echem [-h] [--her E_ads] [--oer OH=..,O=..,OOH=..] [--corrections X=eV] [-U POTENTIAL] [--ph PH] [-T TEMPERATURE] [-o OUTDIR] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono]`

**Optionen:**

| Option | Beschreibung |
|---|---|
| `--her E_ads` | Adsorptionsenergie von H in eV (HER-Reaktion) |
| `--oer OH=..,O=..,OOH=..` | Adsorptionsenergien der drei OER-Intermediate, in eV und bezogen auf Wasser |
| `--corrections X=eV` | thermische Korrekturen ZPE−TΔS pro Intermediat; ohne sie werden die Standardwerte aus der Literatur verwendet |
| `-U, --potential POTENTIAL` | angelegtes Potential in V gegen SHE (bei pH 0 dasselbe wie gegen RHE; der pH rechnet um) |
| `--ph PH` | pH |
| `-T, --temperature TEMPERATURE` | Temperatur in K (Standard: 298.15) |
| `-o, --outdir` | Ausgabeordner (Standard: `echem`) |
| `--dpi DPI` | Auflösung der Bitmap-Abbildungen in Punkten pro Zoll (Standard: 600) |
| `--format` | kommagetrennte Abbildungsformate: pdf,png,svg,eps,tif (Standard: pdf,png) |
| `--no-plot` | nur die Daten exportieren, ohne die Abbildung zu erzeugen |
| `-t, --template` | visuelle Vorlage der Abbildung (die Liste zeigt 'olla-dft templates list') |
| `--font {sans,serif,latex}` | Schriftfamilie: sans, serif oder latex (Computer Modern) |
| `--usetex` | die Texte der Abbildung mit echtem LaTeX setzen (braucht eine LaTeX-Installation) |
| `--palette` | Farbpalette: grayscale, okabe-ito, okabe-ito-dark oder kommagetrennte Hexfarben |
| `--background` | Hintergrundfarbe der Abbildung, z. B. '#FFFFFF' oder 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | Abbildungsbreite je nach Zeitschrift (Standard: generic) |
| `--width` | Breite der Abbildung: single, onehalf oder double, oder eine Zahl in Millimetern |
| `--mono` | Graustufenversion: schwarze Tinte und Linienmuster |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria echem`](THEORY.md)

### `neb`

Reaktionspfad und Aktivierungsbarriere (neb.x)

**Verwendung:** `olla-dft neb [-h] [-o OUTDIR] [--images IMAGES] [--no-ci] [--path-thr PATH_THR] [--nstep NSTEP] [--fix FIX] [--prefix PREFIX] [--collect] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--metal] [--nspin {1,2}] [--mag MAG] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] file [final]`

**Argumente:**

- `file` — Anfangsstruktur (Reaktant)
- `final` — Endstruktur (Produkt)

**Optionen:**

| Option | Beschreibung |
|---|---|
| `-o, --outdir` | Ausgabeordner (Standard: `neb`) |
| `--images IMAGES` | Anzahl der Bilder der Kette (Standard: `7`) |
| `--no-ci` | ohne Climbing Image; die Barriere wird UNTERSCHÄTZT |
| `--path-thr PATH_THR` | Kraftschwelle des Pfads in eV/Å (Standard: `0.05`) |
| `--nstep NSTEP` | maximale Optimierungsschritte des Pfads in neb.x (Standard: 50) |
| `--fix` | Indizes der einzufrierenden Atome (ab 0) |
| `--prefix` | Prefix der Rechnung (wird automatisch erkannt) |
| `--collect` | die Ergebnisse einer bereits gelaufenen Rechnung einlesen, statt die Eingaben vorzubereiten |
| `--pseudo-dir` | Ordner mit den UPF-Pseudopotentialen (ohne Angabe der aus 'olla-dft config') |
| `--pseudo EL=UPF` | ein bestimmtes Pseudopotential erzwingen, z. B. Fe=Fe.rel-pbe.UPF. Kann wiederholt werden. Ohne diese Angabe wählt Olla-DFT mit 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | Cutoff der Wellenfunktionen in Ry (ohne Angabe der von den UPF-Dateien empfohlene) |
| `--ecutrho ECUTRHO` | Cutoff der Dichte in Ry (ohne Angabe der von den UPF-Dateien empfohlene) |
| `--kspacing KSPACING` | Abstand des k-Gitters in Å^-1 |
| `--metal` | metallisches System: Besetzungen mit Smearing statt fester |
| `--nspin {1,2}` | 2 schaltet die Spinpolarisation ein (Standard: 1) |
| `--mag` | Startmagnetisierung: eine Zahl (0.5) oder pro Element (Fe=0.7,O=0). Impliziert --nspin 2 |
| `--dpi DPI` | Auflösung der Bitmap-Abbildungen in Punkten pro Zoll (Standard: 600) |
| `--format` | kommagetrennte Abbildungsformate: pdf,png,svg,eps,tif (Standard: pdf,png) |
| `--no-plot` | nur die Daten exportieren, ohne die Abbildung zu erzeugen |
| `-t, --template` | visuelle Vorlage der Abbildung (die Liste zeigt 'olla-dft templates list') |
| `--font {sans,serif,latex}` | Schriftfamilie: sans, serif oder latex (Computer Modern) |
| `--usetex` | die Texte der Abbildung mit echtem LaTeX setzen (braucht eine LaTeX-Installation) |
| `--palette` | Farbpalette: grayscale, okabe-ito, okabe-ito-dark oder kommagetrennte Hexfarben |
| `--background` | Hintergrundfarbe der Abbildung, z. B. '#FFFFFF' oder 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | Abbildungsbreite je nach Zeitschrift (Standard: generic) |
| `--width` | Breite der Abbildung: single, onehalf oder double, oder eine Zahl in Millimetern |
| `--mono` | Graustufenversion: schwarze Tinte und Linienmuster |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria neb`](THEORY.md)

### `amorphous`

amorpher Festkörper durch Schmelzen und Abschrecken mit einem gelernten Potential

**Verwendung:** `olla-dft amorphous [-h] [-n UNITS] -d G_CM3 [--melt K] [--final K] [--melt-steps MELT_STEPS] [--quench-steps QUENCH_STEPS] [--anneal-steps ANNEAL_STEPS] [--dt FS] [--model MODEL] [--min-dist F] [--seed SEED] [--pack-only] [-o OUTDIR] formula`

**Argumente:**

- `formula` — Formeleinheit, z. B. SiO2

**Optionen:**

| Option | Beschreibung |
|---|---|
| `-n, --units UNITS` | Formeleinheiten in der Zelle (Standard: 8) |
| `-d, --density G_CM3` | Zieldichte in g/cm³ |
| `--melt K` | Schmelztemperatur (Standard: 3000 K) |
| `--final K` | Endtemperatur (Standard: 300 K) |
| `--melt-steps MELT_STEPS` | Dynamikschritte in der Schmelzphase (Standard: 500) |
| `--quench-steps QUENCH_STEPS` | Abschreckschritte: sie legen die Rate fest. Der Standard (1000) ist ein exploratives Abschrecken mit ~3e15 K/s, und der Bericht warnt davor; 27000 senkt sie auf 1e14 K/s |
| `--anneal-steps ANNEAL_STEPS` | Dynamikschritte des Temperns bei der Endtemperatur (Standard: 200) |
| `--dt FS` | Schritt: Integrationsschritt in fs (amorphous) oder Temperaturschritt in K (qha) (Standard: `1.0`) |
| `--model` | interatomares Potential (Standard: `mace`) |
| `--min-dist F` | Faktor auf die Summe der Kovalenzradien beim Packen (Standard: 0.75) |
| `--seed SEED` | Seed; ändern, um eine andere Realisierung zu erzeugen |
| `--pack-only` | nur packen, ohne Dynamik |
| `-o, --outdir` | Ausgabeordner (Standard: `amorfo`) |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria amorphous`](THEORY.md)

## Automatisierung und Qualität

### `doctor`

eine Rechnung diagnostizieren: Konvergenz, Kräfte und warum sie nicht konvergiert

**Verwendung:** `olla-dft doctor [-h] [--system] [--project PROJECT] [--json] [--prefix PREFIX] [-o OUTDIR] [--no-plot] [--dpi DPI] [--format FORMAT] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] [path]`

**Argumente:**

- `path` — Ordner der Rechnung oder Ausgabedatei

**Optionen:**

| Option | Beschreibung |
|---|---|
| `--system` | Installation, Ressourcen, QE und Pseudopotentiale prüfen |
| `--project` | zusätzlich das Quality Gate dieses Projekts prüfen |
| `--json` | die Diagnose als JSON ausgeben |
| `--prefix` | Prefix der Rechnung (wird automatisch erkannt) |
| `-o, --outdir` | Ausgabeordner (Standard: `.`) |
| `--no-plot` | nur die Daten exportieren, ohne die Abbildung zu erzeugen |
| `--dpi DPI` | Auflösung der Bitmap-Abbildungen in Punkten pro Zoll (Standard: 600) |
| `--format` | kommagetrennte Abbildungsformate: pdf,png,svg,eps,tif (Standard: pdf,png) |
| `-t, --template` | visuelle Vorlage der Abbildung (die Liste zeigt 'olla-dft templates list') |
| `--font {sans,serif,latex}` | Schriftfamilie: sans, serif oder latex (Computer Modern) |
| `--usetex` | die Texte der Abbildung mit echtem LaTeX setzen (braucht eine LaTeX-Installation) |
| `--palette` | Farbpalette: grayscale, okabe-ito, okabe-ito-dark oder kommagetrennte Hexfarben |
| `--background` | Hintergrundfarbe der Abbildung, z. B. '#FFFFFF' oder 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | Abbildungsbreite je nach Zeitschrift (Standard: generic) |
| `--width` | Breite der Abbildung: single, onehalf oder double, oder eine Zahl in Millimetern |
| `--mono` | Graustufenversion: schwarze Tinte und Linienmuster |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria doctor`](THEORY.md)

### `audit`

prüfen, dass ein Satz von Rechnungen vergleichbar ist, bevor Energien subtrahiert werden

**Verwendung:** `olla-dft audit [-h] [--index] [--db DB] paths [paths ...]`

**Argumente:**

- `paths` — Ordner oder XML-Dateien der Rechnungen

**Optionen:**

| Option | Beschreibung |
|---|---|
| `--index` | sie zusätzlich in der Datenbank registrieren |
| `--db` | SQLite-Datei des Rechnungsindex (Standard: olla-dft.db) |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria audit`](THEORY.md)

### `crosscheck`

dieselbe Größe über unabhängige Wege kreuzprüfen

**Verwendung:** `olla-dft crosscheck [-h] [-f FILE] [--gap-bandas GAP_BANDAS] [--gap-tauc GAP_TAUC] [project]`

**Argumente:**

- `project` — Projektordner

**Optionen:**

| Option | Beschreibung |
|---|---|
| `-f, --file` | Struktur (für Massen und Volumen) |
| `--gap-bandas GAP_BANDAS` | Bandlücke der Bandstruktur in eV, zur Kreuzprüfung mit der Tauc-Lücke |
| `--gap-tauc GAP_TAUC` | Bandlücke der Tauc-Extrapolation in eV, zur Kreuzprüfung mit der Bandstruktur-Lücke |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria crosscheck`](THEORY.md)

### `cost`

was Olla-DFT über die Geschwindigkeit deiner Maschine weiß

**Verwendung:** `olla-dft cost [-h] [--db DB]`

**Optionen:**

| Option | Beschreibung |
|---|---|
| `--db` | Rechnungsdatenbank (Standard: `olla-dft.db`) |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria cost`](THEORY.md)

### `db`

lokaler Index der Rechnungen

**Verwendung:** `olla-dft db [-h] [--db DB] [-q QUERY] [--export EXPORT] [--formula FORMULA] [--calculation CALCULATION] [--gap-min GAP_MIN] [--gap-max GAP_MAX] [--limit LIMIT] [paths ...]`

**Argumente:**

- `paths` — zu registrierende Ordner

**Optionen:**

| Option | Beschreibung |
|---|---|
| `--db` | SQLite-Datei des Rechnungsindex (Standard: olla-dft.db) |
| `-q, --query` | SQL-Abfrage (nur SELECT) |
| `--export` | alles in eine JSON-Datei exportieren |
| `--formula` | nach Formel filtern, z. B. Si |
| `--calculation` | nach Typ filtern: scf, relax, nscf... |
| `--gap-min GAP_MIN` | minimale Bandlücke in eV |
| `--gap-max GAP_MAX` | maximale Bandlücke in eV |
| `--limit LIMIT` | maximale Anzahl gefilterter Zeilen (Standard: `100`) |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria db`](THEORY.md)

### `hull`

Bildungsenergien und konvexe Hülle

**Verwendung:** `olla-dft hull [-h] [-o OUTDIR] [--elements ELEMENTS] [--threshold THRESHOLD] [--force] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] paths [paths ...]`

**Argumente:**

- `paths` — Ordner oder XML-Dateien der Rechnungen

**Optionen:**

| Option | Beschreibung |
|---|---|
| `-o, --outdir` | Ausgabeordner (Standard: `.`) |
| `--elements` | Reihenfolge der Elemente, z. B. Zn,Al |
| `--threshold THRESHOLD` | Metastabilitätsschwelle in eV/Atom (Standard: `0.025`) |
| `--force` | die Hülle bauen, auch wenn das Audit fehlschlägt |
| `--dpi DPI` | Auflösung der Bitmap-Abbildungen in Punkten pro Zoll (Standard: 600) |
| `--format` | kommagetrennte Abbildungsformate: pdf,png,svg,eps,tif (Standard: pdf,png) |
| `--no-plot` | nur die Daten exportieren, ohne die Abbildung zu erzeugen |
| `-t, --template` | visuelle Vorlage der Abbildung (die Liste zeigt 'olla-dft templates list') |
| `--font {sans,serif,latex}` | Schriftfamilie: sans, serif oder latex (Computer Modern) |
| `--usetex` | die Texte der Abbildung mit echtem LaTeX setzen (braucht eine LaTeX-Installation) |
| `--palette` | Farbpalette: grayscale, okabe-ito, okabe-ito-dark oder kommagetrennte Hexfarben |
| `--background` | Hintergrundfarbe der Abbildung, z. B. '#FFFFFF' oder 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | Abbildungsbreite je nach Zeitschrift (Standard: generic) |
| `--width` | Breite der Abbildung: single, onehalf oder double, oder eine Zahl in Millimetern |
| `--mono` | Graustufenversion: schwarze Tinte und Linienmuster |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria hull`](THEORY.md)

### `mlip`

gelerntes Potential: vorrelaxieren und screenen, bevor DFT-Rechenzeit ausgegeben wird

**Verwendung:** `olla-dft mlip [-h] [-o OUTPUT] [--model {mace,chgnet,m3gnet}] [--size SIZE] [--device DEVICE] [--fmax FMAX] [--steps STEPS] [--fixed-cell] [--span SPAN] [--npoints NPOINTS] [--supercell SUPERCELL] {relax,scan,phonons} file`

**Argumente:**

- `action` {relax,scan,phonons} — auszuführende Aktion (siehe Liste oben)
- `file` — Eingabestruktur (CIF, POSCAR, pw.x-Eingabe...)

**Optionen:**

| Option | Beschreibung |
|---|---|
| `-o, --output` | Ausgabestruktur (relax) |
| `--model {mace,chgnet,m3gnet}` | gelerntes Potential: mace, chgnet oder m3gnet (Standard: mace) |
| `--size` | Größe des MACE-Modells (small/medium/large) (Standard: `small`) |
| `--device` | Gerät, auf dem das Potential läuft: cpu oder cuda (Standard: cpu) |
| `--fmax FMAX` | Zielkraft in eV/Å (Standard: `0.01`) |
| `--steps STEPS` | maximale Relaxationsschritte (Standard: 300) |
| `--fixed-cell` | die Zelle nicht relaxieren, nur die Positionen |
| `--span SPAN` | Bereich des Volumen-Sweeps (scan) (Standard: `0.1`) |
| `--npoints NPOINTS` | Punkte des Volumen-Sweeps (Standard: 15) |
| `--supercell` | Superzelle für das Screening, z. B. 2x2x2 |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria mlip`](THEORY.md)

### `suggest`

Parameter aus deinen früheren Rechnungen vorschlagen

**Verwendung:** `olla-dft suggest [-h] [--db DB] file`

**Argumente:**

- `file` — Eingabestruktur (CIF, POSCAR, pw.x-Eingabe...)

**Optionen:**

| Option | Beschreibung |
|---|---|
| `--db` | SQLite-Datei des Rechnungsindex (Standard: olla-dft.db) |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria suggest`](THEORY.md)

### `datasheet`

Materialdatenblatt und Methodenabsatz

**Verwendung:** `olla-dft datasheet [-h] [-o OUTDIR] [--name NAME] [--methods] [project]`

**Argumente:**

- `project` — Projektordner (Standard: .)

**Optionen:**

| Option | Beschreibung |
|---|---|
| `-o, --outdir` | Ausgabeordner (Standard: `.`) |
| `--name` | Basisname der Dateien |
| `--methods` | nur der Methodenabsatz und die Zitate |

### `report`

lokales Protokoll von Fehlschlägen und Verwirrungen

**Verwendung:** `olla-dft report [-h] [--show SHOW] [--close CLOSE] [--note NOTE] [--stats] [--export EXPORT] [--only-open] [--attach ATTACH] [description ...]`

**Argumente:**

- `description` — was passiert ist (ohne Angabe werden die Vorfälle aufgelistet)

**Optionen:**

| Option | Beschreibung |
|---|---|
| `--show` | einen Vorfall anhand seiner id anzeigen |
| `--close` | einen Vorfall als gelöst markieren |
| `--note` | Notiz beim Schließen |
| `--stats` | welche Unterbefehle am häufigsten fehlschlagen |
| `--export` | alles in eine JSON-Datei packen |
| `--only-open` | nur die offenen Vorfälle auflisten |
| `--attach` | eine Datei anhängen (sie wird ins lokale Protokoll kopiert) |

### `compare`

Läufe vergleichen, ohne inkompatible Energien zu subtrahieren

**Verwendung:** `olla-dft compare [-h] [--reference REFERENCE] [-o OUTPUT] paths [paths ...]`

**Argumente:**

- `paths` — Ordner oder XML-Dateien der Läufe

**Optionen:**

| Option | Beschreibung |
|---|---|
| `--reference REFERENCE` | Index des Referenzlaufs (Standard: 0) |
| `-o, --output` | Vergleich als JSON speichern |

### `tune`

den nächsten Punkt eines Konvergenztests empfehlen

**Verwendung:** `olla-dft tune [-h] [--threshold THRESHOLD] [-o OUTPUT] file`

**Argumente:**

- `file` — CONVERGENCIA.dat

**Optionen:**

| Option | Beschreibung |
|---|---|
| `--threshold THRESHOLD` | Schwelle in meV/Atom (Standard: 1) |
| `-o, --output` | Empfehlung als JSON speichern |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria tune`](THEORY.md)

### `results`

normalisierte Projektergebnisse einlesen, abfragen und exportieren

**Verwendung:** `olla-dft results [-h] [--project PROJECT] [--db DB] [--tag TAG] [--formula FORMULA] [--calculation CALCULATION] [--status {invalid,not_converged,parsed_no_energy,parsed,converged}] [--review-status {unreviewed,accepted,rejected}] [--note NOTE] [--limit LIMIT] [--json] [-o OUTPUT] {ingest,list,show,review,export,explore} [target] [extra_paths ...]`

**Argumente:**

- `action` {ingest,list,show,review,export,explore} — auszuführende Aktion (siehe Liste oben)
- `target` — Eingabepfad für ingest oder id für show
- `extra_paths` — weitere Ordner/XML-Dateien für ingest

**Optionen:**

| Option | Beschreibung |
|---|---|
| `--project` | Projektordner (Standard: .) |
| `--db` | alternative SQLite; standardmäßig .qekit/results.sqlite3 |
| `--tag` | Herkunfts-Tag für ingest |
| `--formula` | nach Formel filtern |
| `--calculation` | nach Rechnungstyp filtern |
| `--status {invalid,not_converged,parsed_no_energy,parsed,converged}` | nach Status filtern: invalid, not_converged, parsed_no_energy, parsed oder converged |
| `--review-status {unreviewed,accepted,rejected}` | bei review Status der menschlichen Prüfung |
| `--note` | bei review Notiz zur Entscheidung |
| `--limit LIMIT` | maximale Datensätze: list=100, explore=10000 |
| `--json` | bei list JSON ausgeben |
| `-o, --output` | Ausgabedatei: export=JSON, explore=interaktives HTML |

### `campaign`

reproduzierbare Matrizen parametrisierter Aufgaben erstellen

**Verwendung:** `olla-dft campaign [-h] [--project PROJECT] [--command CAMPAIGN_COMMAND] [--axis AXIS] [--goal GOAL] [--convergence-file CONVERGENCE_FILE] [--adaptive] [--threshold THRESHOLD] [--execute] [--force] [--parallel PARALLEL] [--retries RETRIES] [--timeout TIMEOUT] [--cancel-file CANCEL_FILE] [-o OUTPUT] {create,list,status,export,run,extend} [target]`

**Argumente:**

- `action` {create,list,status,export,run,extend} — auszuführende Aktion (siehe Liste oben)
- `target` — Name oder id der Kampagne

**Optionen:**

| Option | Beschreibung |
|---|---|
| `--project` | Projektordner (Standard: .) |
| `--command` | Vorlage des olla-dft-Befehls; Felder: {eje}, {index}, {id}, {structure} |
| `--axis` | Achse name=v1,v2; kann wiederholt werden (Standard: `[]`) |
| `--goal` | wissenschaftliches Ziel der Kampagne |
| `--convergence-file` | CONVERGENCIA.dat, aus der eine Empfehlung genommen wird |
| `--adaptive` | den nächsten empfohlenen Wert zur Konvergenzachse hinzufügen |
| `--threshold THRESHOLD` | Konvergenzschwelle beim Erweitern (meV/Atom) |
| `--execute` | bei run die ausgewählten Punkte ausführen |
| `--force` | bei run den Aufgaben-Cache ignorieren |
| `--parallel PARALLEL` | bei run gleichzeitige unabhängige Punkte (Standard: 1) |
| `--retries RETRIES` | Wiederholungen pro fehlgeschlagenem Punkt (Standard: 0) |
| `--timeout TIMEOUT` | maximale Zeit pro Versuch, in Sekunden |
| `--cancel-file` | benutzerdefinierter Marker für kooperativen Abbruch |
| `-o, --output` | JSON-Datei für export |

### `pseudos`

die verfügbaren Pseudopotentiale vergleichen und nach Kriterien wählen, nicht alphabetisch

**Verwendung:** `olla-dft pseudos [-h] [--element ELEMENT] [--task TASK] [--functional FUNCTIONAL] [--cheap] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [file]`

**Argumente:**

- `file` — Struktur

**Optionen:**

| Option | Beschreibung |
|---|---|
| `--element` | kommagetrennte Elemente |
| `--task` | wofür es ist: general, optics, soc, xanes, hubbard, fonones. Jede Aufgabe verwirft die ungeeigneten (Standard: `general`) |
| `--functional` | ein bestimmtes Funktional verlangen (PBE, PZ, PBEsol...) |
| `--cheap` | ultrasoft/PAW bevorzugen, die weniger ebene Wellen brauchen |
| `--pseudo-dir` | Ordner mit den UPF-Pseudopotentialen (ohne Angabe der aus 'olla-dft config') |
| `--pseudo EL=UPF` | ein bestimmtes Pseudopotential erzwingen, z. B. Fe=Fe.rel-pbe.UPF. Kann wiederholt werden. Ohne diese Angabe wählt Olla-DFT mit 'olla-dft pseudos' |

**Wissenschaftliche Grundlagen (Englisch):** [`olla-dft teoria pseudos`](THEORY.md)

## Projekt

### `project`

ein reproduzierbares Projekt verwalten: Quellen, Workflow, Qualität und Dashboard

**Verwendung:** `olla-dft project [-h] [--project PROJECT] [--name NAME] [--command TASK_COMMANDS] [--execute] [--force] [--parallel PARALLEL] [--retries RETRIES] [--timeout TIMEOUT] [--cancel-file CANCEL_FILE] [--reason REASON] [--selftest] [--advanced] [-o OUTPUT] [--pdf] [--theme {auto,light,dark}] [--language {es,en,de}] [--both] [--all-languages] [--verify-environment] [--other OTHER] [--json] {init,add,plan,show,status,validate,run,dashboard,report,export,ingest,environment,diff,cancel,resume} [target]`

**Argumente:**

- `action` {init,add,plan,show,status,validate,run,dashboard,report,export,ingest,environment,diff,cancel,resume} — Aktion auf dem Projekt
- `target` — Verzeichnis, Datei, Ziel, Profil oder Aufgabe je nach Aktion

**Optionen:**

| Option | Beschreibung |
|---|---|
| `--project` | Projekt, von dem aus gearbeitet wird (Standard: .) |
| `--name` | Name beim Initialisieren |
| `--command` | benutzerdefinierte olla-dft-Aufgabe; kann mit plan wiederholt werden |
| `--execute` | run/submit ausführen; standardmäßig wird nur simuliert oder geschrieben |
| `--force` | bei run den Cache ignorieren und alle Aufgaben neu vorbereiten |
| `--parallel PARALLEL` | bei run gleichzeitige unabhängige Aufgaben (Standard: 1) |
| `--retries RETRIES` | Wiederholungen pro fehlgeschlagener Aufgabe (Standard: 0) |
| `--timeout TIMEOUT` | maximale Zeit pro Versuch, in Sekunden |
| `--cancel-file` | benutzerdefinierter Marker für kooperativen Abbruch |
| `--reason` | bei cancel optionaler Grund |
| `--selftest` | bei validate die schnelle Validierung gegen physikalische Referenzen ausführen |
| `--advanced` | bei validate Struktur, Befehle, Einheiten und Kollisionen prüfen |
| `-o, --output` | Ausgabe für dashboard, report oder export |
| `--pdf` | bei report einen eigenständigen PDF-Bericht erzeugen |
| `--theme {auto,light,dark}` | Design des Dashboards (Standard: auto) |
| `--language {es,en,de}` | Oberflächensprache (en, es, de; gemäß --language) |
| `--both` | spanisches und englisches Dashboard in getrennten Dateien erzeugen |
| `--all-languages` | Oberflächen auf Englisch, Spanisch und Deutsch erzeugen |
| `--verify-environment` | bei environment den gespeicherten Lock prüfen |
| `--other` | bei diff Snapshot oder Projekt zum Vergleich |
| `--json` | bei diff JSON ausgeben |

### `resilient`

wiederherstellbare QE-Jobs bei Serverunterbrechungen

**Verwendung:** `olla-dft resilient [-h] [--state STATE] [--pw-cmd PW_CMD] [--runtime-id RUNTIME_ID] [--checkpoint-seconds CHECKPOINT_SECONDS] [--grace-seconds GRACE_SECONDS] [--max-failures MAX_FAILURES] [--threads THREADS] [--keep KEEP] [--max-segments MAX_SEGMENTS] [--resume] [--user USER] [-o OUTPUT] {init,run,status,pause,service} target`

**Argumente:**

- `action` {init,run,status,pause,service} — auszuführende Aktion (siehe Liste oben)
- `target` — Eingabedatei für init; für die übrigen Aktionen das persistente Zustandsverzeichnis des Jobs

**Optionen:**

| Option | Beschreibung |
|---|---|
| `--state` | neues Jobverzeichnis auf einer erhaltenen persistenten Platte |
| `--pw-cmd` | QE- oder MPI-Befehl mit festen Parallelisierungs-Flags (Standard: `pw.x`) |
| `--runtime-id` | Kennung des unveränderlichen Umgebungs-Images |
| `--checkpoint-seconds CHECKPOINT_SECONDS` |  (Standard: `900`) |
| `--grace-seconds GRACE_SECONDS` |  (Standard: `300`) |
| `--max-failures MAX_FAILURES` |  (Standard: `3`) |
| `--threads THREADS` |  (Standard: `1`) |
| `--keep KEEP` | intakte Checkpoint-Generationen, die behalten werden (mindestens 2) (Standard: `2`) |
| `--max-segments MAX_SEGMENTS` | nach dieser Anzahl gespeicherter Segmente anhalten; 0 bedeutet unbegrenzt |
| `--resume` | eine explizite Pause aufheben, bevor fortgesetzt wird |
| `--user` | unprivilegierter Benutzer für den erzeugten systemd-Dienst |
| `-o, --output` | erzeugte Dienstdatei; die Installation erfolgt separat |

## Erscheinungsbild und Konfiguration

### `templates`

Vorlagen auflisten, anzeigen oder exportieren

**Verwendung:** `olla-dft templates [-h] [-o OUTPUT] [{list,show,export}] [name]`

**Argumente:**

- `action` {list,show,export} — list (Standard), show oder export
- `name` — Name der Vorlage

**Optionen:**

| Option | Beschreibung |
|---|---|
| `-o, --output` | JSON-Ausgabedatei (export) |

### `config`

die Konfiguration anzeigen oder ändern

**Verwendung:** `olla-dft config [-h] [{show,set}] [key] [value]`

**Argumente:**

- `action` {show,set} — auszuführende Aktion (siehe Liste oben)
- `key` — Konfigurationsschlüssel, z. B. pseudo_dir, nproc oder language
- `value` — Wert, der dem Schlüssel zugewiesen wird

---

*Olla-DFT 1.5.0*
