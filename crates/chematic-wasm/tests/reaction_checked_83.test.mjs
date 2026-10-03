// Guard checked WASM outcomes and metadata alignment on the exposed 83 rows.
// Full graph/origin/map oracle parity is gated separately by pinned RDKit.
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
  assert.equal(response.products.length, response.product_atom_sources.length, id);
  assert.equal(response.products.length, response.product_template_maps.length, id);
  const reactantSizes = reactants.map((smiles) => {
    try {
      return wasm.parse_smiles(smiles).atom_count();
    } catch {
      return null; // Three jointly invalid inputs include an invalid reactant.
    }
  });
  for (let set = 0; set < response.products.length; set++) {
    assert.equal(response.products[set].length, response.product_atom_sources[set].length, id);
    assert.equal(response.products[set].length, response.product_template_maps[set].length, id);
    for (let product = 0; product < response.products[set].length; product++) {
      const size = wasm.parse_smiles(response.products[set][product]).atom_count();
      const sources = response.product_atom_sources[set][product];
      const maps = response.product_template_maps[set][product];
      assert.equal(sources.length, size, id);
      assert.equal(maps.length, size, id);
      for (let atom = 0; atom < size; atom++) {
        const source = sources[atom];
        if (source !== null) {
          assert.ok(Number.isInteger(source.reactant) && source.reactant >= 0, id);
          assert.ok(Number.isInteger(source.atom) && source.atom >= 0, id);
          assert.ok(source.reactant < reactantSizes.length, id);
          assert.ok(source.atom < reactantSizes[source.reactant], id);
        }
        const map = maps[atom];
        assert.ok(map === null || (Number.isInteger(map) && map > 0 && map <= 65535), id);
      }
    }
  }
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
