import assert from "node:assert/strict";
import { parse_smiles } from "../pkg-node/chematic_wasm.js";

const first = parse_smiles("CC=[C@AL1]=CCC");
const second = parse_smiles("CC=[C@AL2]=CCC");

const firstCanonical = first.canonical_smiles();
const secondCanonical = second.canonical_smiles();
assert.match(firstCanonical, /@AL[12]/);
assert.match(secondCanonical, /@AL[12]/);
assert.notEqual(firstCanonical, secondCanonical);

const firstAgain = parse_smiles(firstCanonical);
const secondAgain = parse_smiles(secondCanonical);
assert.equal(firstAgain.canonical_smiles(), firstCanonical);
assert.equal(secondAgain.canonical_smiles(), secondCanonical);

first.free();
second.free();
firstAgain.free();
secondAgain.free();

console.log("allene axial stereo WASM round-trip: OK");
