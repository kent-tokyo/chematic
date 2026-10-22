import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { createRequire } from 'node:module';
import { fileURLToPath } from 'node:url';
import { createHash } from 'node:crypto';
const require=createRequire(import.meta.url);
const wasm=require(process.argv[2] ?? process.env.CHEMATIC_WASM_MODULE
  ?? fileURLToPath(new URL('../pkg-node/chematic_wasm.js', import.meta.url)));
const packet=JSON.parse(readFileSync(new URL('../../../validation/cjson_mol2_formal_charges/cases.json',import.meta.url)));
for (const c of packet.cases) {
  assert.equal(createHash('sha256').update(c.input).digest('hex'),c.sha256);
  for (const input of [c.input,c.input.replaceAll('\n','\r\n')]) {
    const report=JSON.parse(wasm.convert_with_report_json(input,'cjson','mol2','allow_loss'));
    assert.equal(report.candidate_graph_verified,true,c.id);
    const status=Object.fromEntries(report.observations.map(o=>[o.field,o.status]));
    assert.equal(status.formal_charges,'preserved');
    assert.equal(status.partial_charges,'preserved');
    const after=JSON.parse(wasm.parse_mol2_record_json(report.output));
    assert.deepEqual(after.formal_charges,c.expected.formal_charges);
    assert.deepEqual(after.atoms.map(a=>a.partial_charge),c.expected.partial_charges);
    assert.deepEqual(after.coordinates,c.expected.coords);
    const strict=JSON.parse(wasm.convert_with_report_json(input,'cjson','mol2','strict'));
    assert.equal(strict.output,null);
    assert.equal(strict.candidate_graph_verified,true);
  }
}
const isotope=JSON.parse(packet.cases[0].input);isotope.atoms.isotopes=[18];
const report=JSON.parse(wasm.convert_with_report_json(JSON.stringify(isotope),'cjson','mol2','allow_loss'));
assert.equal(report.candidate_graph_verified,false);
assert.deepEqual(report.candidate_atom_mapping,[]);
assert.equal(report.observations.find(o=>o.field==='isotopes').status,'dropped');
console.log('CJSON/MOL2: 16 indexed formal/partial charge cases, strict refusal and isotope-loss control passed');
