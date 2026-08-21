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
  escaping) — not a naive `split(',')`.
- **SDF** (single or multi-record) — via `chematic-wasm`'s `sdf_to_records_json`. **Hard limits
  inherited from the WASM API, not configurable here**: at most 1,024 records, and the whole
  file capped at 1 MB (whichever limit is hit first — since typical drug-like SDF records run
  2–4 KB, the byte cap usually bites first, around 250–500 records). A record that fails to
  parse is skipped and counted, not silently dropped without mention. This first parsing stage
  is one unavoidable synchronous WASM call (no per-record streaming binding exists yet) — only
  the second stage (computing descriptors per already-parsed record) is chunked/cancelable.
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

- Practical comfort target: a few hundred records processed smoothly.
- A client-side **hard display cap of 2,000 records** applies beyond that — a UI-level safety
  guard, not a WASM API limit. Loading a larger CSV/`.smi` set truncates to the first 2,000 and
  says so.
- Processing (SMILES parsing, descriptor computation, similarity search) is chunked in batches
  of 25 with a yield back to the browser's event loop between chunks, plus a visible progress
  indicator and a Cancel button, so the main thread is never blocked for long on a large input.
- 2D structure thumbnails render lazily (only once a row scrolls into view via
  `IntersectionObserver`), since eagerly depicting hundreds of SVGs on load — not the row count
  itself — is the actual cost driver at this scale.
- No Web Worker is used. This was deliberately not added preemptively; it's a candidate future
  improvement only if main-thread blocking is actually measured to be a problem in practice.

## Known limitations / explicitly out of scope for this version

- No user accounts, cloud save, server-side processing, or `localStorage` auto-save.
- No 3D viewer, structure editor, or reaction editor.
- No machine-learning-based prediction; pKa/ADMET-style scoring is not exposed here at all (the
  Explorer focuses on descriptors/fingerprints/similarity — see the main
  [Playground](../index.html) for pKa/ADMET/3D features).
- Accessibility: keyboard-operable file picker, labelled filter inputs, and pass/fail rendered
  as text (never color alone) are implemented; a full WCAG audit, screen-reader live-region
  progress announcements, and focus-trap modal semantics were **not** attempted in this pass.
- No live-browser verification was performed on this code (no browser automation tool was
  available in the session that built this) — only Node-based unit tests of the pure parsing/
  filtering/sorting/export logic (`tests/explorer.test.mjs`) and a static review of every DOM id
  referenced from `app.js` against the actual markup. A manual pass in a real browser (all input
  paths, filter/sort combinations, similarity search, CSV export/re-import, keyboard-only
  operation, light/dark/mobile) is recommended before relying on this in production.

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
