[**English**](README.md) · [Español](README.es.md)

# Gallery of visual templates

The same silicon bands + DOS figure in every style (`templates`,
`plot -t`).

    olla-dft templates list              # list the available ones
    olla-dft templates show dark         # what each one defines
    olla-dft plot calc/ -t latex-true    # use it

For a template of your own:

    olla-dft templates export dark       # writes an editable JSON
    olla-dft plot calc/ -t dark-copia    # uses it by name

### Files

| File | What it is |
|---|---|
| `galeria_plantillas.png` | side-by-side comparison of all templates |
| `latex.pdf` | Computer Modern without a LaTeX installation |
| `latex-true.pdf` | rendered with real LaTeX (`--usetex`) |
| `dark.pdf` | dark theme for slides |
