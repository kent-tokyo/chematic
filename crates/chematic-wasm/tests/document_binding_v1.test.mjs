import assert from "node:assert/strict";
import path from "node:path";
import { fileURLToPath } from "node:url";

const here = path.dirname(fileURLToPath(import.meta.url));
const root = path.resolve(here, "..", "..", "..");
const wasm = await import(path.join(root, "crates/chematic-wasm/pkg-node/chematic_wasm.js"));

const reaction = {
  id: "rich-rxn-1",
  steps: [{
    id: "step-1",
    components: [
      { id: "reactant-1", role: "reactant", smiles: "[CH3:7]", atom_maps: [{ map_number: 7, atom_index: 0 }], coefficient: 1, origin: "authored" },
      { id: "product-1", role: "product", smiles: "[CH3:7]", atom_maps: [{ map_number: 7, atom_index: 0 }], coefficient: 1, origin: "authored" },
    ],
    conditions: [],
    provenance: [{ source: "fixture", kind: "author", note: "preserve" }],
    origin: "authored",
  }],
  provenance: [{ source: "fixture", kind: "import", note: "preserve" }],
};

const serialized = JSON.parse(wasm.reaction_document_json_v1(JSON.stringify(reaction)));
assert.deepEqual(serialized, reaction);
const editedReaction = JSON.parse(wasm.edit_reaction_document_json_v1(
  JSON.stringify(reaction),
  JSON.stringify({ kind: "set_step_condition", step_id: "step-1", key: "temperature", value: "25 C" }),
));
assert.equal(editedReaction.id, reaction.id);
assert.equal(editedReaction.steps[0].components[0].id, "reactant-1");
assert.deepEqual(editedReaction.provenance, reaction.provenance);
assert.deepEqual(editedReaction.steps[0].provenance, reaction.steps[0].provenance);
assert.deepEqual(editedReaction.steps[0].conditions, [{ key: "temperature", value: "25 C" }]);
assert.throws(() => wasm.reaction_document_json_v1(JSON.stringify({ id: "broken", steps: [] })));
assert.throws(
  () => wasm.reaction_document_to_rxn_v1(JSON.stringify(reaction)),
  /lossy_conversion/,
);

const cdxml = '<CDXML><page id="p1"><arrow id="a1" Custom="keep"/><future id="x"/></page></CDXML>';
const cdxmlEnvelope = JSON.parse(wasm.cdxml_document_json_v1(cdxml));
assert.equal(cdxmlEnvelope.schema, "chematic.cdxml-binding.v1");
assert.equal(cdxmlEnvelope.source, cdxml);
assert.equal(cdxmlEnvelope.document.pages[0].children[0].raw_xml, '<arrow id="a1" Custom="keep"/>');
assert.equal(wasm.cdxml_document_from_json_v1(JSON.stringify(cdxmlEnvelope)), cdxml);
const editedCdxmlEnvelope = JSON.parse(wasm.edit_cdxml_document_json_v1(
  JSON.stringify(cdxmlEnvelope),
  JSON.stringify({ kind: "set_page_attribute", page_id: "p1", key: "title", value: "Page 1" }),
));
const reopened = JSON.parse(wasm.cdxml_document_json_v1(
  wasm.cdxml_document_from_json_v1(JSON.stringify(editedCdxmlEnvelope)),
));
assert.match(reopened.source, /title="Page 1"/);
assert.match(reopened.source, /Custom="keep"/);
assert.equal(reopened.document.pages[0].children[0].kind, "Arrow");
assert.throws(() => wasm.cdxml_document_from_json_v1(JSON.stringify({ schema: "wrong" })));
assert.throws(() => wasm.cdxml_document_projection_json_v1(cdxml), /lossy_conversion/);

console.log("document binding v1 contract: ok");
