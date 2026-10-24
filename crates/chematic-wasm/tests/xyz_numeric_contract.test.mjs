import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { pathToFileURL } from "node:url";
import test from "node:test";

const url = process.env.CHEMATIC_WASM_MODULE
  ? pathToFileURL(process.env.CHEMATIC_WASM_MODULE).href
  : new URL("../pkg-node/chematic_wasm.js", import.meta.url).href;
const wasm = await import(url);
const packet = JSON.parse(readFileSync(new URL("../../../validation/xyz_numeric_read/cases.json", import.meta.url)));
const encode = value => JSON.stringify(value, (_, v) => Object.is(v, -0) ? "__negative_zero__" : v)
  .replaceAll('"__negative_zero__"', "-0.0");
const bits = value => {
  const buffer = new DataView(new ArrayBuffer(8));
  buffer.setFloat64(0, value);
  return buffer.getBigUint64(0).toString(16).padStart(16, "0");
};
function source(token, field) {
  const coords = ["0", "0", "0"];
  if (field < 3) coords[field] = token;
  const header = field === 3 ? `Lattice="1 0 0 0 ${token} 0 0 0 1"`
    : field === 4 ? "Properties=species:S:1:pos:R:3:force:R:3" : "number";
  return `1\n${header}\nC ${coords.join(" ")}${field === 4 ? ` 0 ${token} 0` : ""}\n`;
}

test("finite bits and numeric refusals across topology and JSON readers", () => {
  for (const group of ["valid", "invalid"]) {
    for (const c of packet[group]) {
      for (let field = 0; field < 5; field++) {
        const text = source(c.token, field);
        const readers = field < 3 ? [wasm.mol_from_xyz, wasm.mol_from_extxyz] : [wasm.mol_from_extxyz];
        if (group === "invalid") {
          for (const read of readers) assert.throws(() => read(text), c.id);
          assert.throws(() => wasm.extxyz_frame_json(text), c.id);
        } else {
          for (const read of readers) {
            const mol = read(text);
            assert.equal(mol.atom_count(), 1);
            mol.free();
          }
          const frame = JSON.parse(wasm.extxyz_frame_json(text));
          const value = field < 3 ? frame.coords[0][field]
            : field === 3 ? frame.lattice[4] : frame.properties.force[0][1];
          assert.equal(bits(value), c.bits, c.id);
          if (field === 0) {
            const mol = wasm.mol_from_xyz(text);
            const output = wasm.to_extxyz_json(mol, encode(frame.coords), "{}");
            mol.free();
            assert.equal(bits(JSON.parse(wasm.extxyz_frame_json(output)).coords[0][0]), c.bits, c.id);
          }
        }
      }
    }
  }
});

test("batch JSON records a rejection and retains later finite frames", () => {
  const input = source("-0e-999999", 0) + source("1e-9999", 0) + source("5e-324", 0);
  for (const read of [wasm.xyz_frames_batch_json, wasm.extxyz_frames_batch_json]) {
    const manifest = JSON.parse(read(input, 0, 10));
    assert.equal(manifest.record_count, 3);
    assert.equal(manifest.rejected_count, 1);
    assert.deepEqual(manifest.records.map(r => r.status), ["accepted", "rejected", "accepted"]);
    assert.equal(bits(manifest.records[0].frame.coords[0][0]), "8000000000000000");
    assert.equal(bits(manifest.records[2].frame.coords[0][0]), "0000000000000001");
  }
});

test("string columns and raw metadata retain numeric-looking tokens", () => {
  const text = "1\nProperties=species:S:1:pos:R:3:label:S:1 energy=1e-9999\nC 0 0 0 1e-9999\n";
  const frame = JSON.parse(wasm.extxyz_frame_json(text));
  assert.deepEqual(frame.properties.label, [["1e-9999"]]);
  assert.equal(frame.info.energy, "1e-9999");
});
