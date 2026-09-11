#!/usr/bin/env node

/** Verify the generated Node-WASM artifact's explicit V3000 metadata path. */

import { pathToFileURL } from "node:url";

const packageDir = process.argv[2];
if (!packageDir) throw new Error("usage: node scripts/check_generated_v3000_roundtrip.mjs PKG_DIR");

const wasm = await import(pathToFileURL(`${packageDir}/chematic_wasm.js`));
const block = String.raw`


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
M  V30 BEGIN COLLECTION
M  V30 MDLV30/STEABS ATOMS=(1 1)
M  V30 END COLLECTION
M  V30 BEGIN SGROUP
M  V30 1 SUP 0 ATOMS=(2 1 2)
M  V30 END SGROUP
M  V30 END CTAB
M  END
`;

const roundtripped = wasm.roundtrip_mol_v3000_block(block);
if (!roundtripped.includes("M  V30 1 SUP 0 ATOMS=(2 1 2)")) {
  throw new Error("generated artifact lost V3000 SGROUP metadata");
}
if (roundtripped.indexOf("BEGIN SGROUP") > roundtripped.indexOf("BEGIN COLLECTION")) {
  throw new Error("generated artifact emitted SGROUP after COLLECTION");
}
console.log("generated Node-WASM V3000 roundtrip: OK");
