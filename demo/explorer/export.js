// CSV export for the Local Compound Explorer. csvField/exportToCsv are pure
// and unit-testable; downloadCsv is the one non-pure (Blob/DOM) function.

export const CSV_COLUMNS = [
  'input_index', 'name', 'input_smiles', 'canonical_smiles', 'formula',
  'mw', 'logp', 'tpsa', 'hbd', 'hba', 'rotatable_bonds', 'qed',
  'lipinski_passes', 'pains_passes', 'similarity', 'parse_status', 'error',
];

// Formula-injection guard applies only to string-valued columns whose
// content this app doesn't fully control (user input, or free-text
// error messages) -- never to numeric/boolean columns. A blanket guard
// would break legitimate negative numbers (e.g. logP = -1.03 for
// caffeine), turning a numeric Excel cell into a text cell.
const GUARDED_COLUMNS = new Set([
  'name', 'input_smiles', 'canonical_smiles', 'formula', 'parse_status', 'error',
]);

/**
 * Format one CSV field: applies the OWASP formula-injection guard (a
 * leading apostrophe) to guarded string columns whose value starts with
 * `=`, `+`, `-`, `@`, tab, or CR, then RFC4180-quotes the result if it
 * contains a comma, quote, or newline (independent of the guard above).
 */
export function csvField(value, columnName) {
  let s = value === null || value === undefined ? '' : String(value);
  if (GUARDED_COLUMNS.has(columnName) && /^[=+\-@\t\r]/.test(s)) {
    s = "'" + s;
  }
  if (/[",\n\r]/.test(s)) {
    s = '"' + s.replace(/"/g, '""') + '"';
  }
  return s;
}

function getColumnValue(record, columnName) {
  const d = record.descriptors;
  switch (columnName) {
    case 'input_index': return record.index;
    case 'name': return record.name;
    case 'input_smiles': return record.inputSmiles;
    case 'canonical_smiles': return record.canonicalSmiles ?? '';
    case 'formula': return record.formula ?? '';
    case 'mw': return d ? d.mw : '';
    case 'logp': return d ? d.logP : '';
    case 'tpsa': return d ? d.tpsa : '';
    case 'hbd': return d ? d.hbd : '';
    case 'hba': return d ? d.hba : '';
    case 'rotatable_bonds': return d ? d.rotatableBonds : '';
    case 'qed': return d ? d.qed : '';
    case 'lipinski_passes': return d ? d.lipinskiPasses : '';
    case 'pains_passes': return d ? d.painsPasses : '';
    case 'similarity': return record.similarity ?? '';
    case 'parse_status': return record.status;
    case 'error': return record.errorMessage ?? '';
    default: return '';
  }
}

/** Serialize the given (already filtered+sorted) records as CSV text (CRLF line endings). */
export function exportToCsv(records) {
  const rows = [CSV_COLUMNS.map((c) => csvField(c, c))];
  for (const r of records) {
    rows.push(CSV_COLUMNS.map((c) => csvField(getColumnValue(r, c), c)));
  }
  return rows.map((row) => row.join(',')).join('\r\n') + '\r\n';
}

/** Trigger a browser download of `csvText` as `filename`. Not pure -- DOM/Blob only. */
export function downloadCsv(filename, csvText) {
  const blob = new Blob([csvText], { type: 'text/csv' });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  a.remove();
  URL.revokeObjectURL(url);
}
