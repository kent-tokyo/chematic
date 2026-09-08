import fs from "node:fs";
import * as wasm from "../crates/chematic-wasm/pkg-node/chematic_wasm.js";

const input = fs.readFileSync(0, "utf8");
for (const [index, raw] of input.split(/\r?\n/).entries()) {
  const smiles = raw.trim();
  if (!smiles) continue;
  try {
    const mol = wasm.parse_smiles(smiles);
    console.log(JSON.stringify({ index, smiles, status: "ok", torsion_hex: Buffer.from(wasm.torsion_bitvec(mol)).toString("hex") }));
    mol.free();
  } catch (error) {
    console.log(JSON.stringify({ index, smiles, status: "error", error: String(error) }));
  }
}
