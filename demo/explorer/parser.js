// Pure input-parsing functions for the Local Compound Explorer.
// No DOM access, no WASM calls -- everything here is unit-testable in plain Node.

const SMILES_COLUMN_CANDIDATES = ['smiles', 'SMILES', 'canonical_smiles', 'structure'];
const NAME_COLUMN_CANDIDATES = ['name', 'Name', 'compound', 'id', 'ID'];

/**
 * RFC4180-ish CSV tokenizer: quoted fields may contain commas, embedded
 * newlines, and escaped ("") double quotes. Returns an array of rows, each
 * an array of field strings. Handles CRLF/LF/bare-CR line endings and a
 * missing trailing newline.
 */
export function parseCsvText(text) {
  if (text.charCodeAt(0) === 0xfeff) text = text.slice(1); // strip BOM

  const rows = [];
  let row = [];
  let field = '';
  let state = 'FIELD_START'; // FIELD_START | IN_UNQUOTED | IN_QUOTED | QUOTE_IN_QUOTED

  const endField = () => { row.push(field); field = ''; };
  const endRow = () => { endField(); rows.push(row); row = []; };

  for (let i = 0; i < text.length; i++) {
    const c = text[i];

    if (state === 'QUOTE_IN_QUOTED') {
      if (c === '"') { field += '"'; state = 'IN_QUOTED'; continue; }
      state = 'IN_UNQUOTED'; // the quoted field just closed; reprocess c below
    }

    if (state === 'IN_QUOTED') {
      if (c === '"') { state = 'QUOTE_IN_QUOTED'; }
      else { field += c; }
      continue;
    }

    // FIELD_START or IN_UNQUOTED
    if (c === '"' && state === 'FIELD_START') { state = 'IN_QUOTED'; continue; }
    if (c === ',') { endField(); state = 'FIELD_START'; continue; }
    if (c === '\r') {
      if (text[i + 1] === '\n') i++;
      endRow(); state = 'FIELD_START'; continue;
    }
    if (c === '\n') { endRow(); state = 'FIELD_START'; continue; }
    field += c;
    state = 'IN_UNQUOTED';
  }

  // flush trailing field/row (handles input with no final newline)
  if (field !== '' || row.length > 0 || state !== 'FIELD_START') {
    endRow();
  }

  return rows;
}

/**
 * Case-insensitive exact-match column detection against the fixed
 * candidate lists. Returns {smilesCol, nameCol}, each a 0-based index or
 * null if no header cell matched.
 */
export function detectColumns(headerRow) {
  const lower = headerRow.map((h) => h.trim().toLowerCase());
  const findCol = (candidates) => {
    for (const cand of candidates) {
      const idx = lower.indexOf(cand.toLowerCase());
      if (idx !== -1) return idx;
    }
    return null;
  };
  return {
    smilesCol: findCol(SMILES_COLUMN_CANDIDATES),
    nameCol: findCol(NAME_COLUMN_CANDIDATES),
  };
}

/**
 * Convert parsed CSV rows (header row included) into raw {name, smiles}
 * records, using the given column indices. Skips the header row and any
 * fully blank row. `nameCol` may be null (no name column detected/chosen).
 */
export function csvRowsToRawRecords(rows, smilesCol, nameCol) {
  const records = [];
  for (let i = 1; i < rows.length; i++) {
    const row = rows[i];
    if (row.length === 0 || (row.length === 1 && row[0].trim() === '')) continue;
    const smiles = (row[smilesCol] ?? '').trim();
    if (!smiles) continue;
    const name = nameCol !== null ? (row[nameCol] ?? '').trim() : '';
    records.push({ name, smiles });
  }
  return records;
}

/**
 * Parse `.smi`/newline-SMILES text: one record per non-blank line, SMILES
 * first, an optional whitespace-separated name as the remainder of the
 * line (matches chematic-smiles' own smi_file convention).
 */
export function parseSmiFileText(text) {
  const records = [];
  for (const rawLine of text.split(/\r\n|\r|\n/)) {
    const line = rawLine.trim();
    if (!line) continue;
    const match = line.match(/^(\S+)\s*(.*)$/);
    if (!match) continue;
    records.push({ smiles: match[1], name: match[2].trim() });
  }
  return records;
}

export { SMILES_COLUMN_CANDIDATES, NAME_COLUMN_CANDIDATES };
