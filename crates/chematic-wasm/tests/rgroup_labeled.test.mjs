import assert from "node:assert/strict";
import path from "node:path";
import { fileURLToPath } from "node:url";

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "../../..");
const wasm = await import(path.join(root, "crates/chematic-wasm/pkg-node/chematic_wasm.js"));

const inputs = JSON.stringify(["Cc1ccc(CC)cc1", "CCc1ccc(N)cc1", "CCO"]);
const columns = JSON.parse(
  wasm.rgroup_decompose_columns_json(inputs, "c1cc([*:2])ccc1[*:1]")
);
for (const key of ["Core", "R1", "R2"]) {
  assert.equal(columns[key].length, 3);
  assert.equal(columns[key][2], null);
}
assert.match(columns.Core[0], /\*:1/);
assert.match(columns.Core[0], /\*:2/);
assert.match(columns.R1[0], /\*:1/);
assert.match(columns.R2[0], /\*:2/);

const reordered = JSON.parse(
  wasm.rgroup_decompose_columns_json(
    JSON.stringify(["c1(CC)ccc(C)cc1", "c1(N)ccc(CC)cc1", "CCO"]),
    "c1cc([*:2])ccc1[*:1]"
  )
);
assert.deepEqual(reordered, columns);

assert.throws(
  () => wasm.rgroup_decompose_columns_json(JSON.stringify(["Cc1ccccc1"]), "[c:1]1cccc[c:1]1"),
  /rgroup_duplicate_label:/
);
assert.throws(
  () => wasm.rgroup_decompose_columns_json(JSON.stringify(["Cc1ccccc1"]), "c1cc([*:2]C)ccc1"),
  /rgroup_nonterminal_attachment:/
);

// The historical lowercase row API remains available and unchanged.
const legacy = JSON.parse(wasm.rgroup_decompose_json(JSON.stringify(["Cc1ccccc1"]), "c1ccc(*)cc1"));
assert.equal(legacy[0].matched, true);
assert.equal(typeof legacy[0].r1, "string");

console.log("rgroup_labeled.test.mjs: all assertions passed");
