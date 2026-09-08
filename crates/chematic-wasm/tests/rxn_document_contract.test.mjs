// The Node binding consumes the same typed reaction-document contract as the
// Rust and Python adapters. RXN V2000 itself is intentionally loss-limited.

import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const repoRoot = path.resolve(__dirname, "..", "..", "..");
const wasm = await import(path.join(repoRoot, "crates/chematic-wasm/pkg-node/chematic_wasm.js"));
const fixture = JSON.parse(
  readFileSync(path.join(repoRoot, "validation/cross_binding_contract.json"), "utf8"),
).rxn_document_contract;

assert.equal(fixture.schema_version, 1);
const document = fixture.document;

const rxn = wasm.rxn_document_to_rxn(JSON.stringify(document));
const decoded = JSON.parse(wasm.rxn_document_from_rxn(rxn));
assert.equal(decoded.steps.length, 1);
assert.deepEqual(
  decoded.steps[0].components.map(({ role, smiles }) => ({ role, smiles })),
  fixture.expected_components,
);

assert.throws(
  () => wasm.rxn_document_to_rxn(JSON.stringify({
    ...document,
    steps: [{ ...document.steps[0], conditions: [{ key: "temperature", value: "25 C" }] }],
  })),
);
