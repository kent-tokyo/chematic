import fs from "node:fs";
import * as wasm from "../crates/chematic-wasm/pkg-node/chematic_wasm.js";

const input = fs.readFileSync(0, "utf8");
for (const [index, raw] of input.split(/\r?\n/).entries()) {
  const smiles = raw.trim();
  if (!smiles) continue;
  try {
    const mol = wasm.parse_smiles(smiles);
    const detail = JSON.parse(wasm.rdkit_ecfp4_detail_json(mol));
    const raw_bit_info = Object.entries(detail.rawBitInfo)
      .map(([identifier, entries]) => [Number(identifier), entries.map(([atom, radius]) => [atom, radius]).sort((a, b) => a[0] - b[0] || a[1] - b[1])])
      .sort((a, b) => a[0] - b[0]);
    console.log(JSON.stringify({ index, smiles, status: "ok", raw_bit_info }));
    mol.free();
  } catch (error) {
    console.log(JSON.stringify({ index, smiles, status: "error", error: String(error) }));
  }
}
