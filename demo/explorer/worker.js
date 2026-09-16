// Worker-side molecule analysis for the Local Compound Explorer.
//
// Molecules never cross the structured-clone boundary: this module returns
// plain records only and releases every wasm-bindgen handle before posting a
// result. Rendering remains on the main thread because it creates DOM nodes.

let parseSmiles;
let getDescriptorsJson;
let painsMatchesJson;
let sdfRecordsBatchJson;

async function initialize() {
  const mod = await import("../pkg/chematic_wasm.js");
  const wasmResponse = await fetch("../pkg/chematic_wasm_bg.wasm");
  await mod.default(wasmResponse);
  parseSmiles = mod.parse_smiles;
  getDescriptorsJson = mod.get_descriptors_json;
  painsMatchesJson = mod.pains_matches_json;
  sdfRecordsBatchJson = mod.sdf_records_batch_json;
}

function parseRecord(raw, index) {
  let mol = null;
  const inputIndex = Number.isInteger(raw.index) ? raw.index : index;
  const name = raw.name || `Compound ${inputIndex + 1}`;
  if (raw.sourceError) {
    return {
      index: inputIndex,
      name,
      inputSmiles: raw.smiles || "",
      status: "error",
      canonicalSmiles: null,
      formula: null,
      descriptors: null,
      painsAlerts: [],
      similarity: null,
      errorMessage: raw.sourceError,
    };
  }
  try {
    mol = parseSmiles(raw.smiles);
    return {
      index: inputIndex,
      name,
      inputSmiles: raw.smiles,
      status: "ok",
      canonicalSmiles: mol.canonical_smiles(),
      formula: mol.formula(),
      descriptors: JSON.parse(getDescriptorsJson(mol)),
      painsAlerts: JSON.parse(painsMatchesJson(mol)),
      similarity: null,
      errorMessage: null,
    };
  } catch (error) {
    return {
      index: inputIndex,
      name,
      inputSmiles: raw.smiles,
      status: "error",
      canonicalSmiles: null,
      formula: null,
      descriptors: null,
      painsAlerts: [],
      similarity: null,
      errorMessage: typeof error === "string" ? error : String(error),
    };
  } finally {
    if (mol) {
      try { mol.free(); } catch (_) { /* wasm handle is already invalid */ }
    }
  }
}

self.onmessage = async ({ data }) => {
  try {
    if (data.type === "init") {
      await initialize();
      self.postMessage({ type: "ready", requestId: data.requestId });
      return;
    }
    if (data.type === "parse") {
      const records = data.records.map((raw, offset) => parseRecord(raw, data.startIndex + offset));
      self.postMessage({ type: "parsed", requestId: data.requestId, records });
      return;
    }
    if (data.type === "sdf-batch") {
      const batch = JSON.parse(sdfRecordsBatchJson(data.sdf, data.offset, data.batchSize));
      self.postMessage({ type: "sdf-batch", requestId: data.requestId, batch });
      return;
    }
    throw new Error(`Unsupported explorer worker request: ${data.type}`);
  } catch (error) {
    self.postMessage({
      type: "error",
      requestId: data.requestId,
      message: typeof error === "string" ? error : String(error),
    });
  }
};
