import fs from "node:fs";
import * as wasm from "../crates/chematic-wasm/pkg-node/chematic_wasm.js";

const input = fs.readFileSync(0, "utf8");
for (const [index, raw] of input.split(/\r?\n/).entries()) {
  const smiles = raw.trim();
  if (!smiles) continue;
  try {
    const mol = wasm.parse_smiles(smiles);
    const detail = JSON.parse(wasm.rdkit_ecfp4_detail_json(mol));
    const sparse_counts = Object.entries(detail.sparseCounts)
      .map(([identifier, count]) => [Number(identifier), count])
      .sort((a, b) => a[0] - b[0]);
    console.log(JSON.stringify({ index, smiles, status: "ok", sparse_counts }));
    mol.free();
  } catch (error) {
    console.log(JSON.stringify({ index, smiles, status: "error", error: String(error) }));
  }
}
