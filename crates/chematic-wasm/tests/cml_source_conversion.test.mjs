import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {pathToFileURL} from 'node:url';
const wasm = await import(process.env.CHEMATIC_WASM_MODULE
  ? pathToFileURL(process.env.CHEMATIC_WASM_MODULE).href
  : new URL('../pkg-node/chematic_wasm.js', import.meta.url).href);
const cases = JSON.parse(readFileSync(new URL('../../../validation/cml_source_preservation/cases.json', import.meta.url), 'utf8')).cases;
for (const entry of cases) {
  for (const policy of ['strict', 'allow_loss', 'preserve_source']) {
    const report = JSON.parse(wasm.convert_with_report_json(entry.source, ' .CML ', '.cml', policy));
    assert.equal(report.output, entry.source, entry.id);
    assert.equal(report.status, 'completed');
    assert.equal(report.verification_scope, 'exact_source_bytes');
    assert.equal(report.candidate_graph_verified, false);
    assert.deepEqual(report.candidate_atom_mapping, []);
    assert.deepEqual(report.candidate_bond_mapping, []);
    assert.deepEqual(Object.fromEntries(report.observations.map(o => [o.field, o.status])), {
      source_text: 'preserved', source_semantics: 'opaque_retained',
    });
  }
}
for (const source of ['<molecule/>',
  '<molecule xmlns="http://www.xml-cml.org/schema"><molecule/></molecule>',
  '<molecule xmlns="http://www.xml-cml.org/schema">&#0;</molecule>',
  '<!DOCTYPE molecule><molecule xmlns="http://www.xml-cml.org/schema"/>',
  '<?xml version="1.1"?><molecule xmlns="http://www.xml-cml.org/schema"/>',
  '<molecule xmlns="http://www.xml-cml.org/schema"><atom></molecule>']) {
  assert.throws(() => wasm.convert_with_report_json(source, 'cml', 'cml', 'strict'));
}
assert.throws(() => wasm.convert_with_report_json(cases[0].source, 'cml', 'smiles', 'preserve_source'));
console.log(`CML source conversion: ${cases.length * 3} retained documents and 7 errors passed`);
