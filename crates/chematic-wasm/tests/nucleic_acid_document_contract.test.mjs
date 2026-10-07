import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

const here = path.dirname(fileURLToPath(import.meta.url));
const root = path.resolve(here, "..", "..", "..");
const wasm = await import(path.join(root, "crates/chematic-wasm/pkg-node/chematic_wasm.js"));
const fixture = JSON.parse(
  readFileSync(path.join(root, "validation/nucleic_acid_document_contract.json"), "utf8"),
);

function atPath(value, segments) {
  for (const segment of segments) value = value[segment];
  return value;
}

for (const expected of fixture.valid_cases) {
  const envelope = JSON.parse(wasm.nucleic_acid_validate_json(JSON.stringify(expected.document)));
  assert.equal(envelope.ok, true, expected.id);
  assert.deepEqual(envelope.document, expected.document, expected.id);

  const edited = JSON.parse(
    wasm.nucleic_acid_apply_json_command(
      JSON.stringify(expected.document),
      JSON.stringify(expected.command),
    ),
  );
  assert.equal(edited.ok, true, expected.id);
  assert.deepEqual(atPath(edited.document, expected.edited_path), expected.edited_value, expected.id);
}

for (const expected of fixture.invalid_cases) {
  const envelope = JSON.parse(wasm.nucleic_acid_validate_json(JSON.stringify(expected.document)));
  assert.equal(envelope.ok, false, expected.id);
  assert.equal(envelope.error.code, expected.expected_code, expected.id);
  assert.ok(envelope.error.path.startsWith("/"), expected.id);
}

const overLimit = {
  schema: "chematic.nucleic-acid.v1",
  atom_ids: [],
  strands: [],
  linkages: [],
  annotations: {},
};
for (let index = 0; index < 65; index += 1) {
  const atomId = `a${index}`;
  overLimit.atom_ids.push(atomId);
  overLimit.strands.push({
    id: `s${index}`,
    kind: "dna",
    residues: [{
      id: `r${index}`,
      base: "A",
      sugar: "deoxyribose",
      atom_refs: [atomId],
      annotations: {},
    }],
    annotations: {},
  });
}
const limitEnvelope = JSON.parse(wasm.nucleic_acid_validate_json(JSON.stringify(overLimit)));
assert.equal(limitEnvelope.ok, false);
assert.equal(limitEnvelope.error.code, "resource_limit");
assert.equal(limitEnvelope.error.path, "/strands");
