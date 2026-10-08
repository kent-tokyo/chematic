import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

// The Rust and Python tests replay the same sequences
// (validation/nucleic_acid_edit_sequences.json): every runtime must return
// the same envelopes, a document that validates again after serialization,
// and unchanged atom ownership.
const here = path.dirname(fileURLToPath(import.meta.url));
const root = path.resolve(here, "..", "..", "..");
const wasm = await import(path.join(root, "crates/chematic-wasm/pkg-node/chematic_wasm.js"));
const contract = JSON.parse(
  readFileSync(path.join(root, "validation/nucleic_acid_document_contract.json"), "utf8"),
);
const fixture = JSON.parse(
  readFileSync(path.join(root, "validation/nucleic_acid_edit_sequences.json"), "utf8"),
);

function ownership(document) {
  const owner = new Map();
  for (const strand of document.strands) {
    for (const residue of strand.residues) {
      for (const atom of residue.atom_refs) {
        assert.ok(!owner.has(atom), `atom ${atom} owned twice`);
        owner.set(atom, residue.id);
      }
    }
  }
  for (const atom of document.atom_ids) assert.ok(owner.has(atom), `atom ${atom} unowned`);
  return Object.fromEntries([...owner.entries()].sort());
}

for (const sequence of fixture.sequences) {
  let document = contract.valid_cases.find((c) => c.id === sequence.start).document;
  const owners = ownership(document);
  sequence.steps.forEach((step, k) => {
    const label = `${sequence.id} step ${k}`;
    const envelope = JSON.parse(
      wasm.nucleic_acid_apply_json_command(JSON.stringify(document), JSON.stringify(step.command)),
    );
    assert.equal(envelope.ok, step.ok, label);
    if (step.ok) {
      assert.deepEqual(envelope.document, step.document, label);
      document = envelope.document;
      const reread = JSON.parse(wasm.nucleic_acid_validate_json(JSON.stringify(document)));
      assert.equal(reread.ok, true, label);
      assert.deepEqual(reread.document, document, label);
      assert.deepEqual(ownership(document), owners, label);
    } else {
      assert.equal(envelope.error.code, step.code, label);
      assert.equal(envelope.error.path, step.path, label);
    }
  });
}
