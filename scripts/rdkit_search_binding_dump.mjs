import fs from "node:fs";
import * as wasm from "../crates/chematic-wasm/pkg-node/chematic_wasm.js";

const indexes = new Map();

for (const line of fs.readFileSync(0, "utf8").split(/\r?\n/)) {
  if (!line.trim()) continue;
  const value = JSON.parse(line);
  const key = JSON.stringify(value.db);
  let index = indexes.get(key);
  if (!index) {
    index = new wasm.RdkitSearchIndex(key);
    indexes.set(key, index);
  }
  // The parity runner compares the underlying f64 score. The historical
  // search_json endpoint remains six-decimal truncated for compatibility.
  const output = value.threshold === undefined
    ? index.search_json_precise(value.query, value.k)
    : index.search_json_threshold_precise(value.query, value.threshold, value.k);
  if (output.startsWith("error:")) throw new Error(output);
  console.log(JSON.stringify({ results: JSON.parse(output) }));
}
