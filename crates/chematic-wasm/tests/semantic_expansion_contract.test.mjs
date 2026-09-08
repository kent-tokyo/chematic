import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

const here = path.dirname(fileURLToPath(import.meta.url));
const root = path.resolve(here, "..", "..", "..");
const wasm = await import(path.join(root, "crates/chematic-wasm/pkg-node/chematic_wasm.js"));
const fixture = JSON.parse(
  readFileSync(path.join(root, "validation/cross_binding_contract.json"), "utf8"),
);

assert.equal(fixture.semantic_expansion_contract.schema_version, 1);
for (const expected of fixture.semantic_expansion_contract.cases) {
  const selected = wasm.semantic_apply_json_command(
    JSON.stringify(expected.model),
    JSON.stringify(expected.command),
  );
  const selectedModel = JSON.parse(selected);
  if ("expected_selected_alternative" in expected) {
    assert.equal(
      selectedModel.r_groups[0].selected_alternative,
      expected.expected_selected_alternative,
      expected.id,
    );
  }
  if ("expected_repeat_count" in expected) {
    assert.equal(
      selectedModel.polymer_units[0].repeat_count,
      expected.expected_repeat_count,
      expected.id,
    );
  }
  const expanded = JSON.parse(wasm.semantic_expand_json(expected.base_smiles, selected));
  assert.equal(expanded.schema, "chematic.semantic-expanded.v1", expected.id);
  assert.deepEqual(expanded.source_to_expanded, expected.expected_source_to_expanded, expected.id);
}
