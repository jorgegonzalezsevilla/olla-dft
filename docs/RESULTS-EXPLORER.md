# Result explorer

Version 1.2.0 adds an offline explorer to `project dashboard` and
`results explore`. It reads existing results; it never starts QE or changes an input.

```sh
olla-dft results ingest ./calculation --project ./my-project
olla-dft results explore --project ./my-project -o results.html
# Standalone normalized results database; no project manifest required:
olla-dft --language en results explore --db results.sqlite3 -o results.html
```

1. Open `results.html` in a browser. By default, only converged records are plotted.
2. Filter by state, calculation, formula, label or ID. Select records in the table.
3. Choose numeric X/Y metrics or the record index; expand **Customize presentation**
   for title, color, point/bar style, figure dimensions, text size and Y range.
4. Check the export count. **Export only plotted points** excludes missing metrics,
   different units and values outside the Y range. Unchecking it includes all selected,
   filtered records in CSV, JSON and HTML. SVG/PNG always represent the visible figure.
5. Save SVG for editing, PNG at 1×/2×, CSV for spreadsheets, JSON for exact stored
   values and provenance, or HTML to reopen the selection and presentation offline.

## Scientific interpretation

Values are never converted or rounded in data exports. The table and ticks use
shorter display formatting; hover/focus a point or open its details for stored precision.
When a narrow range requires an axis offset, the figure states the offset explicitly
and labels the delta axis. Record order is not a convergence study or a physical axis.
Bars include zero; use points for a narrow range around a nonzero energy.

Each axis option has an exact unit. A record with another unit is omitted from that
axis, with a visible count. Missing/nonfinite values become JSON null or a blank CSV
cell with a reason. Recorded uncertainty is retained in CSV/JSON; the figure currently
has no error bars. Convergence and human review are separate states.

Composition/method warnings identify exploratory comparisons; matching fingerprints
are not a guarantee of physical comparability. New ingestions include cutoffs (Ry),
k-grid/shift, occupations, spin and smearing. Older records may lack these parameters;
check the original input, pseudopotentials and geometry before interpreting differences.
No automatic deduplication removes distinct result IDs.

## Limits and sharing

The HTML is a snapshot, not a live database connection. Generate it again after new
results arrive. `explore` loads up to 10,000 rows (use `--limit` to reduce it) and displays
loaded/total counts. Graphs require at most 2,000 eligible points; larger selections
must be filtered, with no silent sampling. Tables paginate by 50 rows.

All assets are embedded; there is no external service, font, tracker or browser storage.
An exported HTML embeds only the chosen export scope and its presentation state.
Keep the generated dashboard and its sibling `.results.html` together when copying
that dashboard. Share the explorer alone if you only need the result charts.

Portable exports omit source paths and review notes, but project titles and labels may
contain private text: inspect these before sharing. SHA-256 values identify source
content and recorded methods, not a digital signature. Spreadsheet formula protection
prefixes suspicious text cells with an apostrophe; negative numeric metrics stay numeric.
The CSV uses UTF-8 with BOM, comma delimiters, decimal points and CRLF lines.

This interactive export is separate from the existing `results export` JSON contract,
which remains available for the full normalized records including local provenance.


The record axis follows the snapshot order. CLI/dashboard snapshots record
`ingested_desc_path_asc` in their metadata (latest ingestion first, source path
as a tie-breaker); exported selections retain the original record numbers.
This ordering carries no physical meaning. SVG, JSON and HTML retain the toolkit
version, generation time, schema version and ordering criterion.

For Excel with a decimal-comma locale, import the CSV with **From Text/CSV**,
select comma as the delimiter and a locale that reads the decimal point for metric
columns. Do not rely on double-click detection. In pandas use
`pd.read_csv("olla-results.csv", encoding="utf-8-sig")`. Blank metric cells are
missing values; `<metric>.reason` explains them and `<metric>.unit` supplies the
unit for both value and uncertainty on each row, including mixed-unit exports.
