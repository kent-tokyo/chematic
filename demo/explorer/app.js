import { parseCsvText, detectColumns, csvRowsToRawRecords, parseSmiFileText } from "./parser.js";
import { applyFilters, buildComparator, renderTable } from "./table.js";
import { exportToCsv, downloadCsv } from "./export.js";

const CHUNK_SIZE = 25; // ponytail: fixed constant, tune only if profiling shows it matters
const HARD_RECORD_CAP = 2000; // client-side rendering safety cap, not a WASM API limit

// --- WASM bindings, populated by initWasm() ---
let parseSmiles, getDescriptorsJson, painsMatchesJson, sdfToRecordsJson, tanimotoSmiles, depictSvgOpts, DepictOptions;
let wasmReady = false;

const state = {
  records: [], // CompoundRecord[]
  filters: {},
  sort: { key: "inputOrder", dir: "asc" },
  referenceSmiles: null,
  similarityHasRun: false,
};

let currentAbortController = null;

function $(id) { return document.getElementById(id); }

function showStatus(message) {
  const el = $("explorer-status");
  if (el) el.textContent = message;
}

function showError(message) {
  const el = $("explorer-error");
  if (!el) return;
  el.textContent = message;
  el.classList.remove("hidden");
}

function clearError() {
  const el = $("explorer-error");
  if (el) el.classList.add("hidden");
}

async function initWasm() {
  try {
    const mod = await import("../pkg/chematic_wasm.js");
    const wasmResp = await fetch("../pkg/chematic_wasm_bg.wasm");
    await mod.default(wasmResp);
    parseSmiles = mod.parse_smiles;
    getDescriptorsJson = mod.get_descriptors_json;
    painsMatchesJson = mod.pains_matches_json;
    sdfToRecordsJson = mod.sdf_to_records_json;
    tanimotoSmiles = mod.tanimoto_smiles;
    depictSvgOpts = (mol, opts) => (opts ? mol.depict_svg_opts(opts) : mol.depict_svg());
    DepictOptions = mod.DepictOptions;
    wasmReady = true;
    $("loading-overlay")?.classList.add("hidden");
  } catch (e) {
    showStatus("WASM failed to load: " + String(e));
    throw e;
  }
}

// ---------------------------------------------------------------------------
// Per-record parsing (WASM call + memory-safety discipline)
// ---------------------------------------------------------------------------

function parseOneRecord(raw, index) {
  let mol = null;
  const name = raw.name || `Compound ${index + 1}`;
  try {
    mol = parseSmiles(raw.smiles);
    const canonicalSmiles = mol.canonical_smiles();
    const formula = mol.formula();
    const descriptors = JSON.parse(getDescriptorsJson(mol));
    const painsAlerts = JSON.parse(painsMatchesJson(mol));
    return {
      index, name, inputSmiles: raw.smiles, status: "ok",
      canonicalSmiles, formula, descriptors, painsAlerts,
      similarity: null, errorMessage: null,
    };
  } catch (err) {
    const message = typeof err === "string" ? err : String(err);
    return {
      index, name, inputSmiles: raw.smiles, status: "error",
      canonicalSmiles: null, formula: null, descriptors: null, painsAlerts: [],
      similarity: null, errorMessage: message,
    };
  } finally {
    if (mol) { try { mol.free(); } catch (_) {} }
  }
}

// ---------------------------------------------------------------------------
// Chunked processing driver (progress + cancel, never blocks the main thread)
// ---------------------------------------------------------------------------

async function processRawRecords(rawRecords) {
  clearError();
  if (currentAbortController) currentAbortController.abort();
  const controller = new AbortController();
  currentAbortController = controller;

  const truncated = rawRecords.length > HARD_RECORD_CAP;
  const toProcess = truncated ? rawRecords.slice(0, HARD_RECORD_CAP) : rawRecords;
  if (truncated) {
    showStatus(`Showing the first ${HARD_RECORD_CAP} of ${rawRecords.length} records (client-side display cap).`);
  }

  state.records = [];
  $("explorer-cancel")?.classList.remove("hidden");

  let processed = 0;
  for (let start = 0; start < toProcess.length; start += CHUNK_SIZE) {
    if (controller.signal.aborted) {
      showStatus(`Cancelled after ${processed} of ${toProcess.length} records.`);
      break;
    }
    const chunk = toProcess.slice(start, start + CHUNK_SIZE);
    for (let i = 0; i < chunk.length; i++) {
      state.records.push(parseOneRecord(chunk[i], start + i));
    }
    processed = Math.min(start + CHUNK_SIZE, toProcess.length);
    if (!truncated) showStatus(`Parsing… ${processed}/${toProcess.length}`);
    renderAll();
    await new Promise((resolve) => setTimeout(resolve, 0)); // yield to the event loop
  }

  if (!controller.signal.aborted) {
    const okCount = state.records.filter((r) => r.status === "ok").length;
    const failCount = state.records.length - okCount;
    showStatus(
      failCount === 0
        ? `${okCount} molecule${okCount === 1 ? "" : "s"} loaded.`
        : `${okCount} loaded, ${failCount} failed to parse.`
    );
  }
  $("explorer-cancel")?.classList.add("hidden");
}

// ---------------------------------------------------------------------------
// Similarity search
// ---------------------------------------------------------------------------

async function runSimilaritySearch(referenceSmiles) {
  clearError();
  try {
    const mol = parseSmiles(referenceSmiles);
    mol.free();
  } catch (err) {
    showError("Invalid reference SMILES: " + (typeof err === "string" ? err : String(err)));
    return;
  }

  if (currentAbortController) currentAbortController.abort();
  const controller = new AbortController();
  currentAbortController = controller;

  for (const r of state.records) r.similarity = null;
  state.similarityHasRun = false;
  state.referenceSmiles = referenceSmiles;

  const okRecords = state.records.filter((r) => r.status === "ok");
  $("explorer-cancel")?.classList.remove("hidden");

  for (let start = 0; start < okRecords.length; start += CHUNK_SIZE) {
    if (controller.signal.aborted) break;
    const chunk = okRecords.slice(start, start + CHUNK_SIZE);
    for (const record of chunk) {
      try {
        record.similarity = tanimotoSmiles(referenceSmiles, record.canonicalSmiles);
      } catch (_) {
        record.similarity = null;
      }
    }
    showStatus(`Computing similarity… ${Math.min(start + CHUNK_SIZE, okRecords.length)}/${okRecords.length}`);
    renderAll();
    await new Promise((resolve) => setTimeout(resolve, 0));
  }

  if (!controller.signal.aborted) {
    state.similarityHasRun = true;
    showStatus(`Similarity search complete (ECFP4, radius 2, 2048 bits).`);
    const sortSelect = $("explorer-sort-key");
    if (sortSelect) {
      const opt = sortSelect.querySelector('option[value="similarity"]');
      if (opt) opt.disabled = false;
    }
    const simFilter = $("explorer-filter-similarity-row");
    if (simFilter) simFilter.style.display = "";
  }
  $("explorer-cancel")?.classList.add("hidden");
}

// ---------------------------------------------------------------------------
// Rendering
// ---------------------------------------------------------------------------

function readFiltersFromForm() {
  const val = (id) => {
    const el = $(id);
    if (!el) return "";
    return el.type === "checkbox" ? el.checked : el.value;
  };
  return {
    text: val("explorer-filter-text"),
    validOnly: val("explorer-filter-valid-only"),
    lipinskiPass: val("explorer-filter-lipinski"),
    painsPass: val("explorer-filter-pains"),
    mwMin: val("explorer-filter-mw-min"),
    mwMax: val("explorer-filter-mw-max"),
    logpMin: val("explorer-filter-logp-min"),
    logpMax: val("explorer-filter-logp-max"),
    tpsaMin: val("explorer-filter-tpsa-min"),
    tpsaMax: val("explorer-filter-tpsa-max"),
    qedMin: val("explorer-filter-qed-min"),
    similarityMin: val("explorer-filter-similarity-min"),
  };
}

function renderAll() {
  state.filters = readFiltersFromForm();
  const sortKeyEl = $("explorer-sort-key");
  const sortDirEl = $("explorer-sort-dir");
  if (sortKeyEl) state.sort.key = sortKeyEl.value;
  if (sortDirEl) state.sort.dir = sortDirEl.value;

  const visible = applyFilters(state.records, state.filters)
    .slice()
    .sort(buildComparator(state.sort.key, state.sort.dir));

  const countEl = $("explorer-result-count");
  if (countEl) countEl.textContent = `${visible.length} of ${state.records.length} shown`;

  const emptyEl = $("explorer-empty-filter");
  if (emptyEl) emptyEl.classList.toggle("hidden", visible.length !== 0 || state.records.length === 0);

  const tbody = $("explorer-tbody");
  if (!tbody) return;
  renderTable(
    tbody,
    visible,
    (canonicalSmiles) => {
      let mol = null;
      try {
        mol = parseSmiles(canonicalSmiles);
        const opts = new DepictOptions();
        opts.set_width(160);
        opts.set_height(120);
        const svg = depictSvgOpts(mol, opts);
        opts.free();
        return svg;
      } finally {
        if (mol) mol.free();
      }
    },
    showDetail
  );
}

function showDetail(record) {
  const panel = $("explorer-detail");
  if (!panel) return;
  panel.classList.remove("hidden");

  $("explorer-detail-name").textContent = record.name;
  $("explorer-detail-smiles").textContent = record.canonicalSmiles || record.inputSmiles;

  const svgWrap = $("explorer-detail-svg");
  svgWrap.textContent = "";
  if (record.status === "ok") {
    let mol = null;
    try {
      mol = parseSmiles(record.canonicalSmiles);
      const opts = new DepictOptions();
      opts.set_width(320);
      opts.set_height(240);
      const svgString = depictSvgOpts(mol, opts);
      opts.free();
      const doc = new DOMParser().parseFromString(svgString, "image/svg+xml");
      if (!doc.querySelector("parsererror")) svgWrap.appendChild(doc.documentElement);
    } catch (_) {
      // leave blank on depiction failure
    } finally {
      if (mol) mol.free();
    }
  }

  const descList = $("explorer-detail-descriptors");
  descList.textContent = "";
  if (record.status === "ok") {
    const d = record.descriptors;
    const rows = [
      ["MW", d.mw.toFixed(2)], ["LogP", d.logP.toFixed(2)], ["TPSA", d.tpsa.toFixed(1)],
      ["HBD", d.hbd], ["HBA", d.hba], ["Rotatable bonds", d.rotatableBonds], ["QED", d.qed.toFixed(2)],
      ["Lipinski", d.lipinskiPasses ? "Pass" : "Fail"],
      ["PAINS", d.painsPasses ? "Pass" : `Fail (${record.painsAlerts.join(", ")})`],
      ["Similarity", record.similarity === null ? "—" : record.similarity.toFixed(3) + " (ECFP4, radius 2, 2048 bits)"],
    ];
    for (const [k, v] of rows) {
      const dt = document.createElement("dt"); dt.textContent = k;
      const dd = document.createElement("dd"); dd.textContent = String(v);
      descList.appendChild(dt); descList.appendChild(dd);
    }
  } else {
    const dt = document.createElement("dt"); dt.textContent = "Error";
    const dd = document.createElement("dd"); dd.textContent = record.errorMessage;
    descList.appendChild(dt); descList.appendChild(dd);
  }
}

// ---------------------------------------------------------------------------
// Input handlers
// ---------------------------------------------------------------------------

async function loadSampleDataset() {
  const resp = await fetch("./sample.csv");
  const text = await resp.text();
  const rows = parseCsvText(text);
  const { smilesCol, nameCol } = detectColumns(rows[0]);
  await processRawRecords(csvRowsToRawRecords(rows, smilesCol, nameCol));
}

function loadPastedSmiles(text) {
  const records = parseSmiFileText(text);
  if (records.length === 0) { showError("No SMILES found in the pasted text."); return; }
  processRawRecords(records);
}

function loadCsvText(text) {
  const rows = parseCsvText(text);
  if (rows.length === 0) { showError("Empty CSV file."); return; }
  const { smilesCol, nameCol } = detectColumns(rows[0]);
  if (smilesCol === null) {
    showColumnPicker(rows);
    return;
  }
  processRawRecords(csvRowsToRawRecords(rows, smilesCol, nameCol));
}

function showColumnPicker(rows) {
  const picker = $("explorer-column-picker");
  if (!picker) { showError("Could not detect a SMILES column, and no column picker is available."); return; }
  const select = $("explorer-column-picker-select");
  select.textContent = "";
  rows[0].forEach((header, i) => {
    const opt = document.createElement("option");
    opt.value = String(i);
    opt.textContent = header || `(column ${i + 1})`;
    select.appendChild(opt);
  });
  picker.classList.remove("hidden");
  picker.dataset.pendingRows = "true";
  picker._rows = rows;
}

function loadSdfText(text) {
  let parsed;
  try {
    parsed = JSON.parse(sdfToRecordsJson(text));
  } catch (e) {
    showError("Failed to parse SDF: " + String(e));
    return;
  }
  if (parsed.length === 1 && parsed[0] && parsed[0].error) {
    showError(parsed[0].error);
    return;
  }
  if (parsed.length === 1024) {
    showStatus("SDF reader capped at 1,024 records or ~1 MB input, whichever came first — some records may be missing.");
  }
  const rawRecords = parsed
    .filter((r) => r !== null)
    .map((r) => ({ name: r.name, smiles: r.smiles }));
  const failedCount = parsed.length - rawRecords.length;
  if (failedCount > 0) {
    showStatus(`${failedCount} SDF record(s) failed to parse and were skipped.`);
  }
  processRawRecords(rawRecords);
}

function loadFile(file) {
  const name = file.name.toLowerCase();
  file.text().then((text) => {
    if (name.endsWith(".sdf") || name.endsWith(".mol")) loadSdfText(text);
    else if (name.endsWith(".csv")) loadCsvText(text);
    else if (name.endsWith(".smi") || name.endsWith(".txt")) loadPastedSmiles(text);
    else loadCsvText(text); // best-effort default
  });
}

// ---------------------------------------------------------------------------
// Boot
// ---------------------------------------------------------------------------

function wireEvents() {
  $("explorer-btn-sample")?.addEventListener("click", () => loadSampleDataset());

  $("explorer-btn-parse-paste")?.addEventListener("click", () => {
    loadPastedSmiles($("explorer-paste-textarea").value);
  });

  const fileInput = $("explorer-file-input");
  $("explorer-btn-browse")?.addEventListener("click", () => fileInput.click());
  fileInput?.addEventListener("change", () => {
    if (fileInput.files[0]) loadFile(fileInput.files[0]);
    fileInput.value = "";
  });

  const dropzone = $("explorer-dropzone");
  dropzone?.addEventListener("dragover", (e) => { e.preventDefault(); dropzone.classList.add("drag-over"); });
  dropzone?.addEventListener("dragleave", () => dropzone.classList.remove("drag-over"));
  dropzone?.addEventListener("drop", (e) => {
    e.preventDefault();
    dropzone.classList.remove("drag-over");
    const file = e.dataTransfer.files[0];
    if (file) loadFile(file);
  });

  $("explorer-column-picker-confirm")?.addEventListener("click", () => {
    const picker = $("explorer-column-picker");
    const rows = picker._rows;
    const smilesCol = Number($("explorer-column-picker-select").value);
    picker.classList.add("hidden");
    processRawRecords(csvRowsToRawRecords(rows, smilesCol, null));
  });

  for (const id of [
    "explorer-filter-text", "explorer-filter-valid-only", "explorer-filter-lipinski", "explorer-filter-pains",
    "explorer-filter-mw-min", "explorer-filter-mw-max", "explorer-filter-logp-min", "explorer-filter-logp-max",
    "explorer-filter-tpsa-min", "explorer-filter-tpsa-max", "explorer-filter-qed-min", "explorer-filter-similarity-min",
    "explorer-sort-key", "explorer-sort-dir",
  ]) {
    $(id)?.addEventListener("input", () => renderAll());
    $(id)?.addEventListener("change", () => renderAll());
  }

  $("explorer-btn-reset-filters")?.addEventListener("click", () => {
    document.querySelectorAll("#explorer-filters input[type=text], #explorer-filters input[type=number]")
      .forEach((el) => { el.value = ""; });
    document.querySelectorAll("#explorer-filters input[type=checkbox]").forEach((el) => { el.checked = false; });
    renderAll();
  });

  $("explorer-btn-similarity")?.addEventListener("click", () => {
    const ref = $("explorer-reference-smiles").value.trim();
    if (!ref) { showError("Enter a reference SMILES first."); return; }
    runSimilaritySearch(ref);
  });

  $("explorer-btn-export")?.addEventListener("click", () => {
    const visible = applyFilters(state.records, state.filters)
      .slice()
      .sort(buildComparator(state.sort.key, state.sort.dir));
    downloadCsv("chematic-explorer-export.csv", exportToCsv(visible));
  });

  $("explorer-cancel")?.addEventListener("click", () => {
    if (currentAbortController) currentAbortController.abort();
  });

  $("explorer-detail-close")?.addEventListener("click", () => {
    $("explorer-detail")?.classList.add("hidden");
  });
}

(async () => {
  wireEvents();
  await initWasm();
  showStatus("Ready. Load the sample dataset, paste SMILES, or drop a CSV/SDF/.smi file.");
})();
