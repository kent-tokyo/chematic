import fs from "node:fs";
import * as wasm from "../crates/chematic-wasm/pkg-node/chematic_wasm.js";

const mode = process.argv[2];

function sortedPairs(value) {
  return Object.entries(value)
    .map(([key, entries]) => [
      Number(key),
      entries.map(([atom, radius]) => [atom, radius]).sort((a, b) => a[0] - b[0] || a[1] - b[1]),
    ])
    .sort((a, b) => a[0] - b[0]);
}

const moleculeDumpers = {
  descriptor: (mol) => ({ descriptors: JSON.parse(wasm.get_rdkit_descriptors_json(mol)) }),
  ecfp4: (mol) => ({ ecfp4_hex: Buffer.from(wasm.ecfp4_bitvec(mol)).toString("hex") }),
  maccs: (mol) => ({ maccs_hex: Buffer.from(wasm.maccs_bitvec(mol)).toString("hex") }),
  "rdkit-ecfp4": (mol) => ({ ecfp4_hex: Buffer.from(wasm.rdkit_ecfp4_bitvec(mol)).toString("hex") }),
  "rdkit-ecfp4-bitinfo": (mol) => {
    const detail = JSON.parse(wasm.rdkit_ecfp4_detail_json(mol));
    return { folded_bit_info: sortedPairs(detail.foldedBitInfo) };
  },
  "rdkit-ecfp4-raw-bitinfo": (mol) => {
    const detail = JSON.parse(wasm.rdkit_ecfp4_detail_json(mol));
    return { raw_bit_info: sortedPairs(detail.rawBitInfo) };
  },
  "rdkit-ecfp4-sparse": (mol) => {
    const detail = JSON.parse(wasm.rdkit_ecfp4_detail_json(mol));
    const sparse_counts = Object.entries(detail.sparseCounts)
      .map(([identifier, count]) => [Number(identifier), count])
      .sort((a, b) => a[0] - b[0]);
    return { sparse_counts };
  },
  "rdkit-path": (mol) => ({ path_hex: Buffer.from(wasm.rdkit_path_bitvec(mol)).toString("hex") }),
  "rdkit-rdk": (mol) => ({ rdk_hex: Buffer.from(wasm.rdkit_rdk_bitvec(mol)).toString("hex") }),
  "rdkit-torsion": (mol) => ({ torsion_hex: Buffer.from(wasm.rdkit_torsion_bitvec(mol)).toString("hex") }),
  "topo-path": (mol) => ({ topo_path_hex: Buffer.from(wasm.topo_path_bitvec(mol)).toString("hex") }),
  torsion: (mol) => ({ torsion_hex: Buffer.from(wasm.torsion_bitvec(mol)).toString("hex") }),
};

function dumpSearch(input) {
  const indexes = new Map();
  for (const line of input.split(/\r?\n/)) {
    if (!line.trim()) continue;
    const value = JSON.parse(line);
    const key = JSON.stringify(value.db);
    let index = indexes.get(key);
    if (!index) {
      index = new wasm.RdkitSearchIndex(key);
      indexes.set(key, index);
    }
    const output = value.threshold === undefined
      ? index.search_json_precise(value.query, value.k)
      : index.search_json_threshold_precise(value.query, value.threshold, value.k);
    if (output.startsWith("error:")) throw new Error(output);
    console.log(JSON.stringify({ results: JSON.parse(output) }));
  }
}

function dumpStandardization(input) {
  for (const [index, raw] of input.split(/\r?\n/).entries()) {
    const smiles = raw.trim();
    if (!smiles) continue;
    try {
      const value = wasm.standardize_smiles_report_json(smiles, true, true, true, true);
      if (value.startsWith("error:")) throw new Error(value);
      console.log(JSON.stringify({
        index,
        smiles,
        status: "ok",
        standardized_smiles: JSON.parse(value).smiles,
      }));
    } catch (error) {
      console.log(JSON.stringify({ index, smiles, status: "error", error: String(error) }));
    }
  }
}

function dumpMolecules(input, dumper) {
  for (const [index, raw] of input.split(/\r?\n/).entries()) {
    const smiles = raw.trim();
    if (!smiles) continue;
    let mol;
    try {
      mol = wasm.parse_smiles(smiles);
      console.log(JSON.stringify({ index, smiles, status: "ok", ...dumper(mol) }));
    } catch (error) {
      console.log(JSON.stringify({ index, smiles, status: "error", error: String(error) }));
    } finally {
      mol?.free();
    }
  }
}

const input = fs.readFileSync(0, "utf8");
if (mode === "rdkit-search") {
  dumpSearch(input);
} else if (mode === "standardization") {
  dumpStandardization(input);
} else if (moleculeDumpers[mode]) {
  dumpMolecules(input, moleculeDumpers[mode]);
} else {
  const modes = [...Object.keys(moleculeDumpers), "rdkit-search", "standardization"].join(", ");
  throw new Error(`usage: node scripts/binding_dump.mjs <mode>; modes: ${modes}`);
}
