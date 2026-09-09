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
const application = JSON.parse(
  readFileSync(path.join(repoRoot, "validation/cross_binding_contract.json"), "utf8"),
).reaction_application_contract;
const balance = JSON.parse(
  readFileSync(path.join(repoRoot, "validation/cross_binding_contract.json"), "utf8"),
).reaction_balance_contract;
const center = JSON.parse(
  readFileSync(path.join(repoRoot, "validation/cross_binding_contract.json"), "utf8"),
).reaction_center_contract;

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

const products = JSON.parse(wasm.run_reactants(application.smirks, application.reactants.join("|")));
assert.ok(products.length >= application.expected.minimum_product_sets);
const product = wasm.parse_smiles(products[0][0]);
assert.equal(product.atom_count(), application.expected.product_atom_count);
product.free();
for (const testCase of application.additional_cases) {
  const additionalProducts = JSON.parse(wasm.run_reactants(testCase.smirks, testCase.reactants.join("|")));
  assert.ok(additionalProducts.length >= testCase.minimum_product_sets);
  const additionalProduct = wasm.parse_smiles(additionalProducts[0][0]);
  assert.equal(additionalProduct.atom_count(), testCase.product_atom_count);
  additionalProduct.free();
}
for (const testCase of application.negative_cases) {
  assert.throws(() => wasm.run_reactants(testCase.smirks, testCase.reactants.join("|")));
}

for (const testCase of balance.cases) {
  assert.deepEqual(JSON.parse(wasm.balance_check_json(testCase.reaction)), {
    balanced: testCase.balanced,
    diff: testCase.diff,
  });
}

for (const testCase of center.cases) {
  assert.deepEqual(JSON.parse(wasm.find_reaction_center_json(testCase.reaction)), {
    broken: testCase.broken_bonds,
    formed: testCase.formed_bonds,
    changed: testCase.changed_atoms,
  });
}
