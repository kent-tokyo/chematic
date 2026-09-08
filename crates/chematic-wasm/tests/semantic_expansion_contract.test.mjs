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

const polymerModel = {
  schema: "chematic.semantic.v1",
  atom_ids: ["a1", "a2"],
  bond_ids: [],
  r_groups: [],
  polymer_units: [{
    id: "p1",
    attachment_atoms: ["a1", "a2"],
    end_groups: [],
    repeat_count: null,
    repeat_smiles: "[*]CC[*]",
    repeat_endpoint_atoms: null,
  }],
  extensions: {},
};
const polymerSelected = wasm.semantic_apply_json_command(
  JSON.stringify(polymerModel),
  JSON.stringify({ unit_id: "p1", repeat_count: 3 }),
);
assert.equal(JSON.parse(polymerSelected).polymer_units[0].repeat_count, 3);
const polymerExpanded = JSON.parse(wasm.semantic_expand_json("CC", polymerSelected));
assert.deepEqual(polymerExpanded.source_to_expanded.p1, [2, 3, 4, 5, 6, 7]);
