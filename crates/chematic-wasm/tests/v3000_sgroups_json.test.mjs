import assert from "node:assert/strict";
import path from "node:path";
import { fileURLToPath } from "node:url";

const here = path.dirname(fileURLToPath(import.meta.url));
const root = path.resolve(here, "..", "..", "..");
const wasm = await import(path.join(root, "crates/chematic-wasm/pkg-node/chematic_wasm.js"));

const block = `


  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0 0 0 0
M  V30 2 C 1 0 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 COP 0 ATOMS=(2 1 2) BRKXYZ=(4 1 2 3 4)
M  V30 END SGROUP
M  V30 END CTAB
M  END
`;

const groups = JSON.parse(wasm.v3000_sgroups_json(block));
assert.deepEqual(groups, [{
  id: 1,
  kind: "cop",
  parentId: null,
  atomIds: [1, 2],
  attributes: [{ key: "BRKXYZ", value: "(4 1 2 3 4)" }],
}]);

assert.throws(
  () => wasm.v3000_sgroups_json(block.replace("ATOMS=(2 1 2)", "ATOMS=(2 1)")),
);

console.log("V3000 SGROUP JSON Node/WASM contract: ok");
