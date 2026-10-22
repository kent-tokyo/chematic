import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {pathToFileURL} from 'node:url';
const wasm = await import(process.env.CHEMATIC_WASM_MODULE ? pathToFileURL(process.env.CHEMATIC_WASM_MODULE).href : new URL('../pkg-node/chematic_wasm.js', import.meta.url).href);
const packet = JSON.parse(readFileSync(new URL('../../../validation/gjf-document-2026-10-22/cases.json', import.meta.url), 'utf8'));
for (const c of packet.cases) {
  if (c.expected === 'refuse') {
    assert.throws(() => wasm.roundtrip_gjf_document(c.source), c.id);
    assert.throws(() => wasm.parse_gjf_document_json(c.source), c.id);
    continue;
  }
  assert.deepEqual(Buffer.from(wasm.roundtrip_gjf_document(c.source)), Buffer.from(c.source), c.id);
  const report = JSON.parse(wasm.parse_gjf_document_json(c.source));
  assert.equal(report.schema, 'chematic.gjf-document.v1');
  assert.equal(report.jobs.length, c.jobs.length);
  for (const [index, job] of report.jobs.entries()) {
    assert.equal(job.index, index);
    for (const key of ['source_start','source_end','elements','coords','freeze_codes','charge','multiplicity']) assert.deepEqual(job[key], c.jobs[index][key], `${c.id} ${index} ${key}`);
    assert.ok(!Buffer.from(c.source).subarray(job.source_start,job.source_end).toString().includes('--Link1--'));
  }
}
const one = '# HF/STO-3G\n\nHe\n\n0 1\nHe 0 0 0\n\n';
const exact = Array(128).fill(one).join('--Link1--\n');
assert.equal(wasm.roundtrip_gjf_document(exact), exact);
assert.equal(JSON.parse(wasm.parse_gjf_document_json(exact)).jobs.length,128);
assert.throws(() => wasm.roundtrip_gjf_document(exact+'--Link1--\n'+one));
assert.throws(() => wasm.roundtrip_gjf_record(exact));
assert.throws(() => wasm.convert_common_format(exact,'gjf','com'));
const half = '# HF/STO-3G\n\nHe\n\n0 1\n'+'He 0 0 0\n'.repeat(5000)+'\n';
const atoms = half+'--Link1--\n'+half;
assert.equal(wasm.roundtrip_gjf_document(atoms),atoms);
assert.throws(() => wasm.roundtrip_gjf_document(atoms.trimEnd()+'\nHe 0 0 0\n\n'));
const base = one+'--Link1--\n'+one+'!';
const byteExact = base+'x'.repeat(1_000_000-Buffer.byteLength(base));
assert.equal(wasm.roundtrip_gjf_document(byteExact),byteExact);
assert.throws(() => wasm.roundtrip_gjf_document(byteExact+'x'));
assert.throws(() => wasm.roundtrip_gjf_document(base+'é'.repeat(500000)));
const previous = JSON.parse(readFileSync(new URL('../../../validation/gjf-cartesian-2026-10-22/source-cases.json',import.meta.url),'utf8'));
for (const c of previous.cases) {
  if (c.expected==='retain'||c.id.startsWith('multiple-jobs-')) assert.equal(wasm.roundtrip_gjf_document(c.source),c.source);
  else assert.throws(()=>wasm.roundtrip_gjf_document(c.source));
}
console.log('GJF document: 15 exact/inspected, 22 refused; previous 92 controls, cumulative byte/atom/job boundaries passed');
