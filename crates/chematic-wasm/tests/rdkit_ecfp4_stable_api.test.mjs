// Cross-language shared-corpus test for the promoted RDKit-exact Morgan/ECFP API
// (`rdkit_ecfp4_bitvec` / `rdkit_ecfp4_detail_json` / `rdkit_ecfp_config_bitvec` /
// `rdkit_ecfp_config_detail_json`).
//
// Reads the SAME `validation/ecfp4_rdkit_stable_api_fixtures.json` file the Rust test
// (`crates/chematic-fp/tests/rdkit_morgan_stable_api_fixtures.rs`) and the Python test
// (`crates/chematic-py/tests/test_rdkit_ecfp4_stable_api.py`) read -- generated once
// from a live RDKit oracle by `scripts/gen_ecfp4_rdkit_stable_api_fixtures.py` (RDKit
// version/commit recorded in the file itself).
//
// Not wired into any CI workflow (this repo has no WASM-test CI job at all -- the
// existing `publish-npm.yml`/`pages.yml` only *build* the wasm-pack package). Run
// manually after building the Node target:
//
//   wasm-pack build crates/chematic-wasm --target nodejs --out-dir pkg-node --release
//   node crates/chematic-wasm/tests/rdkit_ecfp4_stable_api.test.mjs

import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import path from "node:path";

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const repoRoot = path.resolve(__dirname, "..", "..", "..");

const wasm = await import(path.join(repoRoot, "crates/chematic-wasm/pkg-node/chematic_wasm.js"));
const { parse_smiles, rdkit_ecfp4_bitvec, rdkit_ecfp4_detail_json, rdkit_ecfp_config_bitvec } = wasm;

const corpus = JSON.parse(
  readFileSync(path.join(repoRoot, "validation/ecfp4_rdkit_stable_api_fixtures.json"), "utf8"),
);

assert.ok(corpus.rdkit_version, "corpus must record the RDKit version");
assert.ok(corpus.fixtures.length >= 30, "expected a real fixture corpus");

function fingerprintBytesToBitList(bytes) {
  const bits = [];
  for (let byteIdx = 0; byteIdx < bytes.length; byteIdx++) {
    const byte = bytes[byteIdx];
    for (let bit = 0; bit < 8; bit++) {
      if (byte & (1 << bit)) bits.push(byteIdx * 8 + bit);
    }
  }
  bits.sort((a, b) => a - b);
  return bits;
}

function sortedPairs(obj) {
  // {key: [[a,r],...]} -> sorted [[key, sorted [[a,r],...]]]
  return Object.entries(obj)
    .map(([k, pairs]) => [Number(k), pairs.map((p) => [p[0], p[1]]).sort((a, b) => a[0] - b[0] || a[1] - b[1])])
    .sort((a, b) => a[0] - b[0]);
}

let checkedOk = 0;
let checkedError = 0;

for (const fx of corpus.fixtures) {
  const mol = parse_smiles(fx.smiles);

  if (fx.expect === "ok") {
    const fpBytes = rdkit_ecfp4_bitvec(mol);
    assert.deepEqual(
      fingerprintBytesToBitList(fpBytes),
      fx.folded_bits,
      `folded_bits mismatch for fixture ${fx.id}`,
    );

    const detail = JSON.parse(rdkit_ecfp4_detail_json(mol));
    assert.deepEqual(fingerprintBytesToBitList(detail.fingerprint), fx.folded_bits);

    const gotCounts = Object.fromEntries(
      Object.entries(detail.sparseCounts).map(([k, v]) => [Number(k), v]),
    );
    const expectedCounts = Object.fromEntries(
      Object.entries(fx.sparse_counts).map(([k, v]) => [Number(k), v]),
    );
    assert.deepEqual(gotCounts, expectedCounts, `sparse_counts mismatch for fixture ${fx.id}`);

    assert.deepEqual(
      sortedPairs(detail.rawBitInfo),
      sortedPairs(fx.raw_bit_info),
      `raw_bit_info mismatch for fixture ${fx.id}`,
    );
    assert.deepEqual(
      sortedPairs(detail.foldedBitInfo),
      sortedPairs(fx.folded_bit_info),
      `folded_bit_info mismatch for fixture ${fx.id}`,
    );

    for (const cell of fx.radius_axis) {
      const bits = fingerprintBytesToBitList(rdkit_ecfp_config_bitvec(mol, cell.radius, 2048));
      assert.deepEqual(bits, cell.folded_bits, `radius axis mismatch for fixture ${fx.id} radius=${cell.radius}`);
    }
    for (const cell of fx.fp_size_axis) {
      const bytes = rdkit_ecfp_config_bitvec(mol, 2, cell.fp_size);
      assert.equal(bytes.length, cell.fp_size / 8);
      const bits = fingerprintBytesToBitList(bytes);
      assert.deepEqual(bits, cell.folded_bits, `fp_size axis mismatch for fixture ${fx.id} fp_size=${cell.fp_size}`);
    }

    checkedOk++;
  } else if (fx.expect === "error") {
    assert.throws(() => rdkit_ecfp4_bitvec(mol), `fixture ${fx.id} expected an error`);
    assert.throws(() => rdkit_ecfp4_detail_json(mol), `fixture ${fx.id} expected an error`);
    checkedError++;
  } else {
    throw new Error(`fixture ${fx.id}: unknown expect '${fx.expect}'`);
  }
}

assert.ok(checkedOk >= 30, `expected real success coverage, got ${checkedOk}`);
assert.ok(checkedError >= 1, `expected real error-path coverage, got ${checkedError}`);

console.log(
  `OK: ${checkedOk} success fixtures + ${checkedError} error fixture(s) matched the shared RDKit oracle corpus ` +
    `(rdkit==${corpus.rdkit_version}) via the WASM bindings.`,
);
