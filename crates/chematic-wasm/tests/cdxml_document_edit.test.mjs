import assert from "node:assert/strict";
import path from "node:path";
import { fileURLToPath } from "node:url";

const here = path.dirname(fileURLToPath(import.meta.url));
const root = path.resolve(here, "..", "..", "..");
const wasm = await import(path.join(root, "crates/chematic-wasm/pkg-node/chematic_wasm.js"));

const input = '<CDXML>\n<page id="p1">\n<arrow id="a1" Custom="keep"/>\n</page>\n<page id="p2">\n<text id="t1"/>\n</page>\n</CDXML>';
const summary = JSON.parse(wasm.cdxml_document_json(input));
assert.equal(summary.pages.length, 2);
assert.equal(summary.pages[1].id, "p2");

const edited = wasm.edit_cdxml_document_json(
  input,
  JSON.stringify({ kind: "set_page_attribute", page_id: "p2", key: "title", value: "Page 2" }),
);
assert.match(edited, /title="Page 2"/);
assert.match(edited, /Custom="keep"/);
