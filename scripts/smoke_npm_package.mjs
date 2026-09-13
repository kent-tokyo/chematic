#!/usr/bin/env node
/** Verify a packed web-target WASM package without relying on Node file fetch. */

import { readFile } from "node:fs/promises";
import { resolve } from "node:path";
import { pathToFileURL } from "node:url";

const packageDir = process.argv[2];
if (!packageDir || process.argv.length !== 3) {
  console.error("usage: node scripts/smoke_npm_package.mjs <installed-package-directory>");
  process.exit(2);
}

const directory = resolve(packageDir);
const moduleUrl = pathToFileURL(resolve(directory, "chematic_wasm.js")).href;
const wasmPath = resolve(directory, "chematic_wasm_bg.wasm");
const { default: init, parse_smiles } = await import(moduleUrl);
const wasmBytes = await readFile(wasmPath);

// The published package is a wasm-pack `web` target. Node's fetch does not
// support file: URLs, so pass the bytes using the generated typed init shape.
await init({ module_or_path: wasmBytes });
const molecule = parse_smiles("c1ccccc1");
const result = { formula: molecule.formula(), atoms: molecule.atom_count() };
molecule.free();

if (result.formula !== "C6H6" || result.atoms !== 6) {
  throw new Error(`unexpected benzene result: ${JSON.stringify(result)}`);
}
console.log(JSON.stringify(result));
