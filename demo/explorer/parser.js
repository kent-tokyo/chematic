// Pure input-parsing functions for the Local Compound Explorer.
// No DOM access, no WASM calls -- everything here is unit-testable in plain Node.

const SMILES_COLUMN_CANDIDATES = ['smiles', 'SMILES', 'canonical_smiles', 'structure'];
const NAME_COLUMN_CANDIDATES = ['name', 'Name', 'compound', 'id', 'ID'];

/**
 * Incremental RFC4180-ish tokenizer. `push()` may end in the middle of an
 * escaped quote, CRLF, or quoted newline; completed rows are returned without
 * retaining preceding input. This lets the Explorer process large CSV files
 * through a `File.stream()` reader rather than first materializing the file.
 */
export class CsvStreamParser {
  constructor() {
    this.row = [];
    this.field = '';
    this.state = 'FIELD_START'; // FIELD_START | IN_UNQUOTED | IN_QUOTED | QUOTE_IN_QUOTED
    this.firstChunk = true;
    this.skipLfAfterCr = false;
  }

  _endField() {
    this.row.push(this.field);
    this.field = '';
  }

  _endRow(rows) {
    this._endField();
    rows.push(this.row);
    this.row = [];
  }

  push(text) {
    const rows = [];
    let start = 0;
    if (this.firstChunk) {
      this.firstChunk = false;
      if (text.charCodeAt(0) === 0xfeff) start = 1; // strip a UTF-8 BOM once
    }
    for (let i = start; i < text.length; i++) {
      const c = text[i];
      if (this.skipLfAfterCr) {
        this.skipLfAfterCr = false;
        if (c === '\n') continue;
      }
      if (this.state === 'QUOTE_IN_QUOTED') {
        if (c === '"') { this.field += '"'; this.state = 'IN_QUOTED'; continue; }
        this.state = 'IN_UNQUOTED'; // reprocess the delimiter/newline below
      }
      if (this.state === 'IN_QUOTED') {
        if (c === '"') this.state = 'QUOTE_IN_QUOTED';
        else this.field += c;
        continue;
      }
      if (c === '"' && this.state === 'FIELD_START') { this.state = 'IN_QUOTED'; continue; }
      if (c === ',') { this._endField(); this.state = 'FIELD_START'; continue; }
      if (c === '\r') {
        this._endRow(rows);
        this.state = 'FIELD_START';
        this.skipLfAfterCr = true;
        continue;
      }
      if (c === '\n') { this._endRow(rows); this.state = 'FIELD_START'; continue; }
      this.field += c;
      this.state = 'IN_UNQUOTED';
    }
    return rows;
  }

  finish() {
    const rows = [];
    // A closing quote is valid at EOF; an unclosed quote is retained as field
    // content, matching the historical best-effort parser rather than guessing
    // a different row boundary.
    if (this.field !== '' || this.row.length > 0 || this.state !== 'FIELD_START') {
      this._endRow(rows);
    }
    this.state = 'FIELD_START';
    this.skipLfAfterCr = false;
    return rows;
  }
}

/** Parse a complete CSV string through the same incremental tokenizer. */
export function parseCsvText(text) {
  const parser = new CsvStreamParser();
  return [...parser.push(text), ...parser.finish()];
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
