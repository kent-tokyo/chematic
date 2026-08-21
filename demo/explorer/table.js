// Record-list filtering, sorting, and DOM rendering for the Local Compound
// Explorer. applyFilters/buildComparator/matchesFreeText are pure and
// unit-testable; renderTable is the only DOM-touching function here.

/**
 * True if `record` matches a case-insensitive free-text query against its
 * name, input SMILES, canonical SMILES, or formula. An empty/blank query
 * always matches.
 */
export function matchesFreeText(record, query) {
  const q = query.trim().toLowerCase();
  if (!q) return true;
  const haystacks = [record.name, record.inputSmiles, record.canonicalSmiles, record.formula];
  return haystacks.some((h) => typeof h === 'string' && h.toLowerCase().includes(q));
}

function inRange(value, min, max) {
  if (value === null || value === undefined) return false;
  if (min !== null && min !== undefined && min !== '' && value < Number(min)) return false;
  if (max !== null && max !== undefined && max !== '' && value > Number(max)) return false;
  return true;
}

/**
 * Apply the combined filter set to a record array. Every filter is
 * independently optional (null/undefined/'' means "not active").
 */
export function applyFilters(records, filters) {
  const f = filters || {};
  return records.filter((r) => {
    if (f.validOnly && r.status !== 'ok') return false;
    if (!matchesFreeText(r, f.text || '')) return false;

    if (r.status === 'ok') {
      const d = r.descriptors;
      if (f.lipinskiPass && !d.lipinskiPasses) return false;
      if (f.painsPass && !d.painsPasses) return false;
      if (!inRangeIfSet(d.mw, f.mwMin, f.mwMax)) return false;
      if (!inRangeIfSet(d.logP, f.logpMin, f.logpMax)) return false;
      if (!inRangeIfSet(d.tpsa, f.tpsaMin, f.tpsaMax)) return false;
      if (f.qedMin !== null && f.qedMin !== undefined && f.qedMin !== '' && d.qed < Number(f.qedMin)) {
        return false;
      }
    } else if (f.lipinskiPass || f.painsPass || hasAnyNumericFilter(f)) {
      // an error row can never satisfy a descriptor-based filter
      return false;
    }

    if (f.similarityMin !== null && f.similarityMin !== undefined && f.similarityMin !== '') {
      if (r.similarity === null || r.similarity === undefined) return false;
      if (r.similarity < Number(f.similarityMin)) return false;
    }

    return true;
  });
}

function inRangeIfSet(value, min, max) {
  const minSet = min !== null && min !== undefined && min !== '';
  const maxSet = max !== null && max !== undefined && max !== '';
  if (!minSet && !maxSet) return true;
  return inRange(value, minSet ? min : null, maxSet ? max : null);
}

function hasAnyNumericFilter(f) {
  return ['mwMin', 'mwMax', 'logpMin', 'logpMax', 'tpsaMin', 'tpsaMax', 'qedMin'].some(
    (k) => f[k] !== null && f[k] !== undefined && f[k] !== ''
  );
}

const SORT_ACCESSORS = {
  name: (r) => (r.name || '').toLowerCase(),
  mw: (r) => (r.status === 'ok' ? r.descriptors.mw : null),
  logP: (r) => (r.status === 'ok' ? r.descriptors.logP : null),
  tpsa: (r) => (r.status === 'ok' ? r.descriptors.tpsa : null),
  qed: (r) => (r.status === 'ok' ? r.descriptors.qed : null),
  similarity: (r) => r.similarity,
  inputOrder: (r) => r.index,
};

/**
 * Build a comparator for `records.sort(...)`. `sortKey` is one of the
 * SORT_ACCESSORS keys; `sortDir` is 'asc' or 'desc'. Records whose sort
 * value is null/undefined always sort to the end, regardless of
 * direction (matters for `similarity` before any search has run, and for
 * numeric columns on error rows).
 */
export function buildComparator(sortKey, sortDir) {
  const accessor = SORT_ACCESSORS[sortKey] || SORT_ACCESSORS.inputOrder;
  const dir = sortDir === 'desc' ? -1 : 1;
  return (a, b) => {
    const va = accessor(a);
    const vb = accessor(b);
    const aNull = va === null || va === undefined;
    const bNull = vb === null || vb === undefined;
    if (aNull && bNull) return 0;
    if (aNull) return 1;
    if (bNull) return -1;
    if (va < vb) return -dir;
    if (va > vb) return dir;
    return 0;
  };
}

function svgToElement(svgString) {
  const parser = new DOMParser();
  const doc = parser.parseFromString(svgString, 'image/svg+xml');
  if (doc.querySelector('parsererror')) return null;
  return doc.documentElement;
}

/**
 * Render the visible rows into `tbody`. `depictFn(canonicalSmiles)` is
 * called lazily (only when a row's SVG cell scrolls into view) and must
 * return raw SVG markup; it is inserted via DOMParser, never innerHTML.
 */
export function renderTable(tbody, visibleRecords, depictFn, onRowSelect) {
  tbody.textContent = '';
  const fragment = document.createDocumentFragment();

  for (const record of visibleRecords) {
    const tr = document.createElement('tr');
    tr.dataset.index = String(record.index);

    const svgCell = document.createElement('td');
    svgCell.className = 'explorer-svg-cell';
    svgCell.setAttribute('aria-label', record.status === 'ok' ? `Structure of ${record.name}` : 'Parse error');
    tr.appendChild(svgCell);

    const nameCell = document.createElement('td');
    nameCell.textContent = record.name;
    tr.appendChild(nameCell);

    if (record.status === 'ok') {
      const d = record.descriptors;
      for (const val of [
        record.formula,
        d.mw.toFixed(2),
        d.logP.toFixed(2),
        d.tpsa.toFixed(1),
        d.hbd,
        d.hba,
        d.rotatableBonds,
        d.qed.toFixed(2),
        d.lipinskiPasses ? 'Pass' : 'Fail',
        d.painsPasses ? 'Pass' : 'Fail',
        record.similarity === null || record.similarity === undefined ? '—' : record.similarity.toFixed(3),
      ]) {
        const td = document.createElement('td');
        td.textContent = String(val);
        tr.appendChild(td);
      }
    } else {
      const td = document.createElement('td');
      td.colSpan = 11; // matches the 11 data cells the 'ok' branch above appends
      td.className = 'explorer-error-cell';
      td.textContent = record.errorMessage || 'Parse error';
      tr.appendChild(td);
    }

    tr.tabIndex = 0;
    tr.addEventListener('click', () => onRowSelect(record));
    tr.addEventListener('keydown', (e) => {
      if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); onRowSelect(record); }
    });

    fragment.appendChild(tr);
  }

  tbody.appendChild(fragment);

  if (typeof IntersectionObserver === 'undefined') return;
  const rows = tbody.querySelectorAll('tr');
  const observer = new IntersectionObserver((entries, obs) => {
    for (const entry of entries) {
      if (!entry.isIntersecting) continue;
      const tr = entry.target;
      obs.unobserve(tr);
      const index = Number(tr.dataset.index);
      const record = visibleRecords.find((r) => r.index === index);
      if (!record || record.status !== 'ok') continue;
      const svgCell = tr.querySelector('.explorer-svg-cell');
      if (!svgCell) continue;
      try {
        const svgEl = svgToElement(depictFn(record.canonicalSmiles));
        if (svgEl) svgCell.appendChild(svgEl);
      } catch (_) {
        // leave the cell empty on depiction failure -- not fatal to the row
      }
    }
  }, { root: null, rootMargin: '200px' });
  rows.forEach((tr) => observer.observe(tr));
}
