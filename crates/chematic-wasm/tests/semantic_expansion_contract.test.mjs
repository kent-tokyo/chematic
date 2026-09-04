import assert from "node:assert/strict";
import path from "node:path";
import { fileURLToPath } from "node:url";

const here = path.dirname(fileURLToPath(import.meta.url));
const root = path.resolve(here, "..", "..", "..");
const wasm = await import(path.join(root, "crates/chematic-wasm/pkg-node/chematic_wasm.js"));

const model = {
  schema: "chematic.semantic.v1",
  atom_ids: ["a1", "a2"],
  bond_ids: [],
  r_groups: [{
    id: "r1",
    attachment_atoms: ["a2"],
    alternatives: ["[*]O"],
    selected_alternative: null,
  }],
  polymer_units: [],
  extensions: {},
};

const selected = wasm.semantic_apply_json_command(
  JSON.stringify(model),
  JSON.stringify({ group_id: "r1", alternative: 0 }),
);
const expanded = JSON.parse(wasm.semantic_expand_json("CC", selected));
assert.equal(expanded.schema, "chematic.semantic-expanded.v1");
assert.deepEqual(expanded.source_to_expanded.r1, [2]);
