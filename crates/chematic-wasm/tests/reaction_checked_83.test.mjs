// Guard the checked WASM outcome contract on the exposed 83-row reaction set.
// Graph, atom-origin and template-map parity are checked by the Rust source
// audit, not by this JSON-only adapter test.
import assert from "node:assert/strict";
import { createHash } from "node:crypto";
import { readFileSync } from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "../../..");
const wasm = await import(path.join(root, "crates/chematic-wasm/pkg-node/chematic_wasm.js"));
const readFixture = (relativePath) => {
  const bytes = readFileSync(path.join(root, relativePath));
  return {
    data: JSON.parse(bytes.toString("utf8")),
    sha256: createHash("sha256").update(bytes).digest("hex"),
  };
};

const base = readFixture("validation/reaction_product_parity_cases.json");
const strata = readFixture("validation/reaction_product_parity_strata_v2.json");
const audit = readFixture("validation/results/v1.0.31-source-reaction-83-checked-summary.json").data;
assert.equal(base.sha256, audit.fixtures.base_sha256);
assert.equal(strata.sha256, audit.fixtures.strata_sha256);
assert.equal(strata.data.base_fixture_sha256, base.sha256);
assert.equal(audit.compatibility_profile, "rdkit-2026.03.6");

const cases = [...base.data, ...strata.data.cases];
assert.equal(cases.length, 83);
assert.equal(audit.accounting.input, cases.length);
const expected = new Map(audit.rows.map(({ id, outcome }) => [id, outcome]));
assert.equal(expected.size, cases.length);

const counts = new Map();
for (const { id, smirks, reactants } of cases) {
  const outcome = expected.get(id);
  assert.ok(outcome, `missing audited outcome: ${id}`);
  const response = JSON.parse(wasm.run_reactants_checked(smirks, reactants.join("|"), true));
  assert.equal(response.profile, "rdkit-2026.03.6", id);
  assert.ok(Array.isArray(response.products), id);
  counts.set(response.status, (counts.get(response.status) ?? 0) + 1);
  switch (outcome) {
    case "graph_origin_map_match":
      assert.ok(["products", "no_match"].includes(response.status), id);
      assert.equal(response.reason, null, id);
      break;
    case "typed_unsupported":
      assert.equal(response.status, "typed_unsupported", id);
      assert.equal(response.reason, "chiral_reactant_template_semantics", id);
      assert.deepEqual(response.products, [], id);
      break;
    case "typed_or_diagnosed_refusal":
      assert.equal(response.status, "typed_refusal", id);
      assert.equal(response.reason, "product_valence", id);
      assert.ok(response.valence_rejected_matches > 0, id);
      break;
    case "joint_invalid_input":
      assert.equal(response.status, "typed_refusal", id);
      assert.ok(["smiles_parse", "smirks_parse", "reactant_count_mismatch"].includes(response.reason), id);
      break;
    default:
      assert.fail(`unknown audited outcome for ${id}: ${outcome}`);
  }
}
assert.deepEqual(Object.fromEntries(counts), {
  products: 69,
  no_match: 7,
  typed_unsupported: 3,
  typed_refusal: 4,
});
console.log("checked WASM reaction outcomes: 83 exposed rows classified");
