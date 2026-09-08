import fs from "node:fs";
import * as wasm from "../crates/chematic-wasm/pkg-node/chematic_wasm.js";

const input = fs.readFileSync(0, "utf8");
for (const [index, raw] of input.split(/\r?\n/).entries()) {
  const smiles = raw.trim();
  if (!smiles) continue;
  try {
    const value = wasm.standardize_smiles_report_json(smiles, true, true, true, true);
    if (value.startsWith("error:")) throw new Error(value);
    const report = JSON.parse(value);
    console.log(JSON.stringify({ index, smiles, status: "ok", standardized_smiles: report.smiles }));
  } catch (error) {
    console.log(JSON.stringify({ index, smiles, status: "error", error: String(error) }));
  }
}
