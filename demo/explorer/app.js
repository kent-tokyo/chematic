import { CsvStreamParser, parseCsvText, detectColumns, csvRowsToRawRecords, parseSmiFileText } from "./parser.js";
import { applyFilters, buildComparator, renderTable } from "./table.js";
import { exportToCsv, downloadCsv } from "./export.js";
import { formatBatchOutcome, summarizeKnownLengthBatch } from "./batch-accounting.js";

// 100 records keeps the Worker message overhead bounded at the workflow
// cap while still yielding to the event loop between batches. Individual
// records are small, bounded WASM operations; this is not a hidden unbounded
// batch size.
const CHUNK_SIZE = 100;
const SDF_BATCH_SIZE = 64;
// Keep the complete analysed record set available for filtering, search, and
// export, while MAX_RENDERED_ROWS bounds the DOM. This is a product workflow
// bound, separate from the much smaller per-call WASM API bounds.
const HARD_RECORD_CAP = 100_000;
const SMALL_INPUT_RENDER_PROGRESS_STRIDE = 500;
const LARGE_INPUT_RENDER_PROGRESS_STRIDE = 5_000;
const MAX_RENDERED_ROWS = 250;

// --- WASM bindings, populated by initWasm() ---
let parseSmiles, getDescriptorsJson, painsMatchesJson, tanimotoSmiles, depictSvgOpts, DepictOptions;
let wasmReady = false;
let analysisWorker = null;
let nextWorkerRequestId = 1;
const pendingWorkerRequests = new Map();
// A Worker request cannot be preempted once posted.  This monotonically
// increasing generation prevents its late response from changing a newer
// import after cancellation or replacement.
let activeImportGeneration = 0;

const state = {
  records: [], // CompoundRecord[]
  filters: {},
  sort: { key: "inputOrder", dir: "asc" },
  referenceSmiles: null,
  similarityHasRun: false,
  // Last known-length import outcome.  Incomplete work retains its unprocessed
  // range instead of manufacturing skipped records or a successful summary.
  batchOutcome: null,
};

let currentAbortController = null;

function $(id) { return document.getElementById(id); }

function renderProgressStride(totalRecords) {
  return totalRecords > 10_000
    ? LARGE_INPUT_RENDER_PROGRESS_STRIDE
    : SMALL_INPUT_RENDER_PROGRESS_STRIDE;
}

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
    tanimotoSmiles = mod.tanimoto_smiles;
    depictSvgOpts = (mol, opts) => (opts ? mol.depict_svg_opts(opts) : mol.depict_svg());
    DepictOptions = mod.DepictOptions;
    wasmReady = true;
    $("loading-overlay")?.classList.add("hidden");
  } catch (e) {
    // The loading overlay is a full-screen, high-z-index element that would
    // otherwise permanently hide this message behind "Loading chematic
    // (WASM)…" -- surface the failure where it's actually visible.
    const overlay = $("loading-overlay");
    if (overlay) overlay.textContent = "Failed to load chematic (WASM): " + String(e);
    showStatus("WASM failed to load: " + String(e));
    throw e;
  }
}

// ---------------------------------------------------------------------------
// Worker transport. The main thread owns the display-only WASM module while
// parsing and descriptor work runs in a separately initialized Worker.
// ---------------------------------------------------------------------------

function workerRequest(type, payload = {}) {
  const requestId = nextWorkerRequestId++;
  return new Promise((resolve, reject) => {
    pendingWorkerRequests.set(requestId, { resolve, reject });
    analysisWorker.postMessage({ type, requestId, ...payload });
  });
}

function closeAnalysisWorker() {
  if (!analysisWorker) return;
  const error = new Error("Explorer Worker closed.");
  for (const { reject } of pendingWorkerRequests.values()) reject(error);
  pendingWorkerRequests.clear();
  analysisWorker.terminate();
  analysisWorker = null;
}

async function initAnalysisWorker() {
  if (typeof Worker === "undefined") throw new Error("Web Worker is unavailable in this browser.");
  analysisWorker = new Worker(new URL("./worker.js", import.meta.url), { type: "module" });
  analysisWorker.onmessage = ({ data }) => {
    const pending = pendingWorkerRequests.get(data.requestId);
    if (!pending) return;
    pendingWorkerRequests.delete(data.requestId);
    if (data.type === "error") pending.reject(new Error(data.message));
    else pending.resolve(data);
  };
  analysisWorker.onerror = (event) => {
    const error = new Error(event.message || "Explorer Worker failed.");
    for (const { reject } of pendingWorkerRequests.values()) reject(error);
    pendingWorkerRequests.clear();
  };
  await workerRequest("init");
  document.documentElement.dataset.explorerAnalysis = "worker";
}

async function parseChunkInWorker(records, startIndex) {
  const response = await workerRequest("parse", { records, startIndex });
  return response.records;
}

async function readSdfBatchInWorker(sdf, offset) {
  const response = await workerRequest("sdf-batch", { sdf, offset, batchSize: SDF_BATCH_SIZE });
  return response.batch;
}

// ---------------------------------------------------------------------------
// Chunked processing driver (progress + cancel, never blocks the main thread)
// ---------------------------------------------------------------------------

async function processRawRecords(rawRecords) {
  clearError();
  if (currentAbortController) currentAbortController.abort();
  const controller = new AbortController();
  currentAbortController = controller;
  const importGeneration = ++activeImportGeneration;

  const truncated = rawRecords.length > HARD_RECORD_CAP;
  const toProcess = truncated ? rawRecords.slice(0, HARD_RECORD_CAP) : rawRecords;
  const progressStride = renderProgressStride(toProcess.length);
  if (truncated) {
    showStatus(`Showing the first ${HARD_RECORD_CAP} of ${rawRecords.length} records (client-side display cap).`);
  }

  state.records = [];
  state.batchOutcome = null;
  $("explorer-cancel")?.classList.remove("hidden");

  let processed = 0;
  let terminalReason = truncated ? "client_record_cap" : null;
  for (let start = 0; start < toProcess.length; start += CHUNK_SIZE) {
    if (controller.signal.aborted) {
      terminalReason = "cancelled";
      break;
    }
    const chunk = toProcess.slice(start, start + CHUNK_SIZE);
    const records = await parseChunkInWorker(chunk, start);
    if (controller.signal.aborted || importGeneration !== activeImportGeneration) {
      // Cancellation can arrive while a Worker request is in flight.  Report
      // the same terminal state as the pre-request cancellation path instead
      // of silently falling out of the loop.
      terminalReason = "cancelled";
      break;
    }
    state.records.push(...records);
    processed = Math.min(start + CHUNK_SIZE, toProcess.length);
    if (!truncated) showStatus(`Parsing… ${processed}/${toProcess.length}`);
    // Preserve all completed records for filtering/search/export, but avoid
    // rebuilding thousands of table rows for every worker reply.
    if (processed % progressStride === 0 || processed === toProcess.length) renderAll();
    await new Promise((resolve) => setTimeout(resolve, 0)); // yield to the event loop
  }

  // A later import owns the UI state. Its caller has already initialized a
  // new accounting result, so an older response must not overwrite records,
  // status text, or the cancel control.
  if (importGeneration !== activeImportGeneration) return;

  const outcome = summarizeKnownLengthBatch({
    inputCount: rawRecords.length,
    completedCount: state.records.length,
    terminalReason: processed === rawRecords.length ? null : terminalReason || "cancelled",
  });
  state.batchOutcome = outcome;
  const okCount = state.records.filter((r) => r.status === "ok").length;
  showStatus(formatBatchOutcome(outcome, {
    acceptedCount: okCount,
    rejectedCount: state.records.length - okCount,
  }));
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

// Re-run parse/descriptor work for failed SMILES without discarding their
// original order, names, or source text. This is useful after a transient
// Worker failure; malformed SDF source records intentionally remain evidence
// rows and are not offered for retry because they have no SMILES payload.
async function retryFailedSmiles() {
  const retryable = state.records.filter(
    (record) => record.status === "error" && Boolean(record.inputSmiles?.trim())
  );
  if (retryable.length === 0) {
    showStatus("No failed SMILES rows can be retried. SDF source rejections are preserved as evidence.");
    return;
  }

  clearError();
  if (currentAbortController) currentAbortController.abort();
  const controller = new AbortController();
  currentAbortController = controller;
  const positions = new Map(state.records.map((record, position) => [record.index, position]));
  let completed = 0;
  let recovered = 0;

  state.referenceSmiles = null;
  state.similarityHasRun = false;
  const similarityOption = $("explorer-sort-key")?.querySelector('option[value="similarity"]');
  if (similarityOption) similarityOption.disabled = true;
  const similarityFilter = $("explorer-filter-similarity-row");
  if (similarityFilter) similarityFilter.style.display = "none";
  $("explorer-cancel")?.classList.remove("hidden");

  try {
    for (let start = 0; start < retryable.length; start += CHUNK_SIZE) {
      if (controller.signal.aborted) break;
      const chunk = retryable.slice(start, start + CHUNK_SIZE).map((record) => ({
        index: record.index,
        name: record.name,
        smiles: record.inputSmiles,
      }));
      const replacements = await parseChunkInWorker(chunk, 0);
      if (controller.signal.aborted) break;
      for (const replacement of replacements) {
        const position = positions.get(replacement.index);
        if (position === undefined) continue;
        if (replacement.status === "ok") recovered += 1;
        state.records[position] = replacement;
      }
      completed += replacements.length;
      showStatus(`Retrying failed SMILES… ${completed}/${retryable.length}`);
      renderAll();
      await new Promise((resolve) => setTimeout(resolve, 0));
    }
    if (!controller.signal.aborted) {
      const stillFailed = retryable.length - recovered;
      showStatus(`Retry complete: ${recovered} recovered; ${stillFailed} still failed.`);
    } else {
      showStatus(`Retry cancelled after ${completed}/${retryable.length} failed SMILES rows.`);
    }
  } finally {
    if (controller === currentAbortController) $("explorer-cancel")?.classList.add("hidden");
  }
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
  if (countEl) {
    const rendered = Math.min(visible.length, MAX_RENDERED_ROWS);
    countEl.textContent = rendered === visible.length
      ? `${visible.length} of ${state.records.length} shown`
      : `${rendered} rendered of ${visible.length} matching (${state.records.length} loaded)`;
  }

  const retryButton = $("explorer-btn-retry-failed");
  if (retryButton) {
    const retryable = state.records.some(
      (record) => record.status === "error" && Boolean(record.inputSmiles?.trim())
    );
    retryButton.classList.toggle("hidden", !retryable);
  }

  const emptyEl = $("explorer-empty-filter");
  if (emptyEl) emptyEl.classList.toggle("hidden", visible.length !== 0 || state.records.length === 0);

  const tbody = $("explorer-tbody");
  if (!tbody) return;
  renderTable(
    tbody,
    visible.slice(0, MAX_RENDERED_ROWS),
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

async function loadCsvFile(file) {
  // The normal auto-detected CSV path is streaming: rows are sent to the Worker
  // in bounded chunks and the whole input is never copied into one JS string.
  // A file with an unknown SMILES column deliberately falls back to the existing
  // manual-picker path, where retaining rows is necessary for the user to choose.
  if (typeof file.stream !== "function") {
    // Older browsers keep the historical, bounded-by-workflow fallback instead
    // of failing an otherwise valid local CSV import.
    loadCsvText(await file.text());
    return;
  }
  const reader = file.stream().getReader();
  const decoder = new TextDecoder();
  const csv = new CsvStreamParser();
  let header = null;
  let smilesCol = null;
  let nameCol = null;
  let pending = [];
  let processed = 0;
  let sawDataRow = false;
  let truncated = false;
  let controller = null;

  const start = () => {
    if (currentAbortController) currentAbortController.abort();
    controller = new AbortController();
    currentAbortController = controller;
    state.records = [];
    clearError();
    $("explorer-cancel")?.classList.remove("hidden");
  };
  const flush = async () => {
    if (!controller || controller.signal.aborted || pending.length === 0) return;
    const chunk = pending;
    pending = [];
    const records = await parseChunkInWorker(chunk, processed);
    if (controller.signal.aborted) return;
    state.records.push(...records);
    processed += records.length;
    showStatus(`Reading CSV… ${processed} molecule${processed === 1 ? "" : "s"} analysed.`);
    if (processed % renderProgressStride(processed) === 0) renderAll();
    await new Promise((resolve) => setTimeout(resolve, 0));
  };
  const consumeRows = async (rows) => {
    for (const row of rows) {
      if (header === null) {
        header = row;
        ({ smilesCol, nameCol } = detectColumns(header));
        if (smilesCol === null) return false;
        start();
        continue;
      }
      if (row.length === 0 || (row.length === 1 && row[0].trim() === "")) continue;
      sawDataRow = true;
      const smiles = (row[smilesCol] ?? "").trim();
      if (!smiles) continue;
      if (processed + pending.length >= HARD_RECORD_CAP) {
        truncated = true;
        return true;
      }
      pending.push({ name: nameCol !== null ? (row[nameCol] ?? "").trim() : "", smiles });
      if (pending.length >= CHUNK_SIZE) await flush();
      if (controller?.signal.aborted) return true;
    }
    return true;
  };

  try {
    for (;;) {
      const { value, done } = await reader.read();
      const text = value ? decoder.decode(value, { stream: !done }) : "";
      const accepted = await consumeRows(csv.push(text));
      if (header !== null && smilesCol === null) {
        await reader.cancel();
        // Preserve the manual-column-picker behavior only when auto-detection
        // cannot decide the schema. It is intentionally not the large-file path.
        loadCsvText(await file.text());
        return;
      }
      if (!accepted || truncated || controller?.signal.aborted || done) break;
    }
    if (!truncated && !controller?.signal.aborted) await consumeRows(csv.finish());
    if (header === null || !sawDataRow) {
      showError("Empty CSV file.");
      return;
    }
    if (controller?.signal.aborted) {
      showStatus(`Cancelled after ${processed} CSV records.`);
      return;
    }
    await flush();
    renderAll();
    const failures = state.records.filter((record) => record.status !== "ok").length;
    const base = failures === 0
      ? `${processed} molecule${processed === 1 ? "" : "s"} loaded.`
      : `${processed - failures} loaded, ${failures} failed to parse.`;
    showStatus(truncated ? `${base} Showing the first ${HARD_RECORD_CAP} records (client-side display cap).` : base);
  } catch (error) {
    if (!controller?.signal.aborted) {
      state.records = [];
      renderAll();
      showError("Failed to read CSV: " + (typeof error === "string" ? error : String(error)));
    }
  } finally {
    if (controller === currentAbortController) $("explorer-cancel")?.classList.add("hidden");
  }
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

async function loadSdfText(text) {
  clearError();
  if (currentAbortController) currentAbortController.abort();
  const controller = new AbortController();
  currentAbortController = controller;
  state.records = [];
  $("explorer-cancel")?.classList.remove("hidden");

  let offset = 0;
  let accepted = 0;
  let rejected = 0;
  try {
    for (;;) {
      if (controller.signal.aborted) {
        showStatus(`Cancelled after SDF record ${offset}. ${accepted} loaded, ${rejected} rejected.`);
        break;
      }
      const batch = await readSdfBatchInWorker(text, offset);
      if (!batch || !Array.isArray(batch.records) || !Number.isInteger(batch.next_offset)) {
        throw new Error("SDF batch returned an invalid manifest.");
      }
      if (batch.offset !== offset || batch.next_offset < offset) {
        throw new Error("SDF batch did not advance deterministically.");
      }
      const rawRecords = batch.records.map((entry) => {
        const inputIndex = entry.input_index;
        if (entry.status === "accepted" && entry.record) {
          return { index: inputIndex, name: entry.record.name, smiles: entry.record.smiles };
        }
        rejected += 1;
        return {
          index: inputIndex,
          name: `SDF record ${inputIndex + 1}`,
          smiles: "",
          sourceError: "SDF record was rejected before descriptor analysis.",
        };
      });
      const records = await parseChunkInWorker(rawRecords, offset);
      if (controller.signal.aborted) {
        showStatus(`Cancelled after SDF record ${offset}. ${accepted} loaded, ${rejected} rejected.`);
        break;
      }
      state.records.push(...records);
      accepted += records.filter((record) => record.status === "ok").length;
      offset = batch.next_offset;
      showStatus(`Reading SDF… ${offset} records inspected; ${accepted} loaded, ${rejected} rejected.`);
      if (offset % renderProgressStride(offset) === 0 || batch.status === "complete") renderAll();
      if (batch.status === "complete") break;
      if (batch.status !== "partial" || batch.next_offset === batch.offset) {
        throw new Error("SDF batch did not provide a resumable continuation.");
      }
      await new Promise((resolve) => setTimeout(resolve, 0));
    }
    if (!controller.signal.aborted) {
      showStatus(`${accepted} loaded, ${rejected} rejected from ${offset} SDF records.`);
    }
  } catch (error) {
    state.records = [];
    renderAll();
    showError("Failed to read SDF: " + (typeof error === "string" ? error : String(error)));
  } finally {
    $("explorer-cancel")?.classList.add("hidden");
  }
}

async function loadFile(file) {
  const name = file.name.toLowerCase();
  if (name.endsWith(".csv")) {
    await loadCsvFile(file);
    return;
  }
  const text = await file.text();
  if (name.endsWith(".sdf") || name.endsWith(".mol")) void loadSdfText(text);
  else if (name.endsWith(".smi") || name.endsWith(".txt")) loadPastedSmiles(text);
  else loadCsvText(text); // best-effort default
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
    if (fileInput.files[0]) void loadFile(fileInput.files[0]);
    fileInput.value = "";
  });

  const dropzone = $("explorer-dropzone");
  dropzone?.addEventListener("dragover", (e) => { e.preventDefault(); dropzone.classList.add("drag-over"); });
  dropzone?.addEventListener("dragleave", () => dropzone.classList.remove("drag-over"));
  dropzone?.addEventListener("drop", (e) => {
    e.preventDefault();
    dropzone.classList.remove("drag-over");
    const file = e.dataTransfer.files[0];
    if (file) void loadFile(file);
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

  $("explorer-btn-retry-failed")?.addEventListener("click", () => {
    void retryFailedSmiles();
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
  await Promise.all([initWasm(), initAnalysisWorker()]);
  showStatus("Ready. Load the sample dataset, paste SMILES, or drop a CSV/SDF/.smi file.");
})();

window.addEventListener("pagehide", closeAnalysisWorker, { once: true });
