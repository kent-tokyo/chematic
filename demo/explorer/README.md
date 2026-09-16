# Local Compound Explorer

A static, browser-only tool to load, view, filter, sort, similarity-search, and CSV-export a
batch of compounds — entirely client-side, powered by the `chematic-wasm` build already used by
the [Playground](../index.html). No server, no build step beyond the existing WASM bundle.

Live: <https://kent-tokyo.github.io/chematic/explorer/>

## Supported input formats

- **CSV** — auto-detects a SMILES column (`smiles`, `SMILES`, `canonical_smiles`, `structure`)
  and a name column (`name`, `Name`, `compound`, `id`, `ID`) by exact, case-insensitive header
  match. If no SMILES column is found, a manual column picker is shown. The parser is a
  hand-written RFC4180-style tokenizer (quoted fields, embedded commas/newlines, doubled-quote
  escaping) — not a naive `split(',')`. Auto-detected files stream from `File.stream()` into
  bounded Worker chunks; the rare manual-column-picker fallback reads the file so it can retain
  rows while the user selects a column.
- **SDF** (single or multi-record) — via the resumable `sdf_records_batch_json` WASM API in
  fixed 64-record Worker batches. **Hard limits inherited from the WASM API, not configurable
  here**: the whole file is capped at 1 MB and each request at 1,024 records. An invalid SDF
  record remains in the result set as an error row with its original input index; it is neither
  silently dropped nor counted as a successfully analysed molecule. Cancelling stops before the
  next batch, leaving completed rows available for filtering and export; loading the file again
  starts the deterministic sequence from record zero.
- **`.smi` / newline-separated SMILES / pasted text** — one record per non-blank line, SMILES
  first, an optional whitespace-separated name as the rest of the line.
- **Built-in sample dataset** (`sample.csv`) — 16 well-known drugs/small molecules, verified
  against the real WASM build (correct formula + molecular weight for every entry).

## Computed descriptors

Per molecule: canonical SMILES, molecular formula, MW, LogP (Crippen), TPSA, HBD, HBA,
rotatable bonds, QED, Lipinski pass/fail, PAINS pass/fail (plus which PAINS alerts fired) — all
from `chematic-wasm`'s existing `get_descriptors_json`/`pains_matches_json` exports, not
recomputed in JavaScript.

## Similarity search

Uses `tanimoto_smiles` (ECFP4, radius 2, 2048 bits — chematic's own fingerprint hash, not
RDKit-bit-identical) against a reference SMILES you provide. This routes similarity search
entirely through SMILES strings with no `MolHandle` ever created for it, which is why it needs
no separate memory-management discipline from the main parse path.

## CSV export

Exports the *current filtered + sorted view* (not just the loaded set) with columns:
`input_index, name, input_smiles, canonical_smiles, formula, mw, logp, tpsa, hbd, hba,
rotatable_bonds, qed, lipinski_passes, pains_passes, similarity, parse_status, error`.

**Formula-injection protection**: string-valued columns (`name`, `input_smiles`,
`canonical_smiles`, `formula`, `parse_status`, `error`) that start with `=`, `+`, `-`, `@`, tab,
or CR are prefixed with a leading apostrophe (the standard OWASP CSV-injection guard), which
forces spreadsheet software to treat the cell as text rather than evaluating it as a formula.
This guard is deliberately **not** applied to numeric/boolean columns (`mw`, `logp`, `tpsa`,
etc.) — a legitimate negative value like LogP = −1.03 must stay a real number, not get turned
into a text cell.

## Privacy model

Supported analysis runs entirely in your browser's own WASM sandbox. Uploaded/pasted molecule
data is never sent to a server by this page. Reloading the page clears the session — there is
no auto-save to `localStorage` in this version. This page makes no network calls of its own
beyond loading its own static assets (WASM binary, CSS, this HTML/JS).

## Performance / limits

- Parsing and descriptor analysis run in a dedicated ES-module **Web Worker**.
  Only DOM rendering and on-demand depiction stay on the main thread. The worker
  imports the same-origin WASM asset explicitly, so CSP must allow `worker-src
  'self'` and `script-src 'self'`; initialization failure is surfaced instead
  of silently falling back to synchronous main-thread processing. Pending work
  is rejected and the Worker is terminated on `pagehide`.
- Practical comfort target: a few hundred records processed smoothly. The Trust
  Release gate is a separate 10,000-record Worker workflow; this demo's current
  rendering cap remains intentionally lower. CSV and SMI are the 10,000-record
  workflow inputs; the bounded 1 MB SDF API is intentionally a separate contract.
- A client-side **workflow cap of 10,000 records** applies beyond that — a UI-level safety
  guard, not a WASM API limit. Loading a larger CSV/`.smi` set truncates to the first 10,000 and
  says so. All accepted records remain available to filtering, similarity and CSV export; the
  table deliberately renders at most 250 at a time to keep interaction responsive.
- Worker parsing/descriptor computation and main-thread similarity search are processed in
  batches of 100 with a yield back to the browser event loop between chunks, plus a visible
  progress indicator and a Cancel button.
- Auto-detected CSV input uses the same 100-record Worker batches while incrementally decoding
  the file. An RFC4180 quoted field, escaped quote, or CRLF split across byte chunks remains one
  logical field/row; the parser does not assume chunk boundaries are line boundaries.
- CSV export includes accepted and rejected rows. It follows the active filter and sort; the
  default `Input order` sort preserves the original 0-based input index.
- Failed rows that contain SMILES can be retried without replacing the dataset; the retry keeps
  their input order, name, source text, and export row. Source-level SDF rejections have no
  SMILES payload and remain visible evidence rather than being retried.
- 2D structure thumbnails render lazily (only once a row scrolls into view via
  `IntersectionObserver`), since eagerly depicting hundreds of SVGs on load — not the row count
  itself — is the actual cost driver at this scale.

## Known limitations / explicitly out of scope for this version

- No user accounts, cloud save, server-side processing, or `localStorage` auto-save.
- No 3D viewer, structure editor, or reaction editor.
- No machine-learning-based prediction; pKa/ADMET-style scoring is not exposed here at all (the
  Explorer focuses on descriptors/fingerprints/similarity — see the main
  [Playground](../index.html) for pKa/ADMET/3D features).
- Accessibility: keyboard-operable file picker, labelled filter inputs, and pass/fail rendered
  as text (never color alone) are implemented; a full WCAG audit, screen-reader live-region
  progress announcements, and focus-trap modal semantics were **not** attempted in this pass.
- Browser smoke coverage includes `scripts/explorer_worker_10k_smoke.mjs` and
  `scripts/explorer_worker_10k_smoke.py` for the 10,000-row SMI path. The Python smoke also
  covers the equivalent CSV-file path, checks Cancel, and disables network access after
  page/Worker initialization before exercising those local paths. It downloads the full
  10,000-row CSV export and verifies input order/status rather than accepting the 250-row
  render window as evidence. `scripts/explorer_sdf_worker_smoke.py` covers resumable SDF batches with an
  inline rejected record and verifies its exported input order/status. `scripts/explorer_retry_failed_smoke.py`
  verifies that a malformed-SMILES retry retains the rejected export row. All three scripts accept
  `--engine chromium|firefox|webkit`; Chromium additionally requires `--browser /path/to/browser`.
  The Python scripts require a local static server.
  Manual checks of filters, similarity, CSV re-import, keyboard-only use, light/dark themes, and
  mobile layout remain necessary before relying on this in production.

## Files

```
demo/explorer/
  index.html    DOM skeleton only
  app.js        WASM loading, chunked processing, event wiring, DOM rendering
  parser.js     Pure: CSV tokenizer, column detection, .smi parsing
  table.js      Pure: filter/sort/free-text logic, plus the one DOM-rendering function
  export.js     Pure: CSV field escaping + serialization, plus the one Blob/download function
  styles.css
  sample.csv
  package.json  {"type": "module"} only, so Node resolves the test file's ES imports correctly
  tests/explorer.test.mjs   node:assert unit tests for every pure function above
```
