// The Node binding consumes the same typed reaction-document contract as the
// Rust and Python adapters. RXN V2000 itself is intentionally loss-limited.

import assert from "node:assert/strict";
import path from "node:path";
import { fileURLToPath } from "node:url";

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const repoRoot = path.resolve(__dirname, "..", "..", "..");
const wasm = await import(path.join(repoRoot, "crates/chematic-wasm/pkg-node/chematic_wasm.js"));

const document = {
  id: "rxn-contract",
  steps: [{
    id: "step-1",
    components: [
      { id: "reactant-1", role: "reactant", smiles: "CC", coefficient: 1, origin: "authored" },
      { id: "product-1", role: "product", smiles: "CC", coefficient: 1, origin: "authored" },
    ],
    conditions: [],
    provenance: [],
    origin: "authored",
  }],
  provenance: [],
};

const rxn = wasm.rxn_document_to_rxn(JSON.stringify(document));
const decoded = JSON.parse(wasm.rxn_document_from_rxn(rxn));
assert.equal(decoded.steps.length, 1);
assert.deepEqual(
  decoded.steps[0].components.map(({ role, smiles }) => ({ role, smiles })),
  [
    { role: "reactant", smiles: "CC" },
    { role: "product", smiles: "CC" },
  ],
);

assert.throws(
  () => wasm.rxn_document_to_rxn(JSON.stringify({
    ...document,
    steps: [{ ...document.steps[0], conditions: [{ key: "temperature", value: "25 C" }] }],
  })),
);
