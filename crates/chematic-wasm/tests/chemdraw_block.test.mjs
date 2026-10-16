import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {pathToFileURL} from 'node:url';
const wasm=await import(process.env.CHEMATIC_WASM_MODULE?pathToFileURL(process.env.CHEMATIC_WASM_MODULE).href:new URL('../pkg-node/chematic_wasm.js',import.meta.url).href);
const fixture=JSON.parse(readFileSync(new URL('../../../validation/rdkit_issues/chemdraw-block/cases.json',import.meta.url),'utf8'));
for (const c of fixture.cases) {
  const data=Buffer.from(c.cdx_hex,'hex');
  for (const mode of ['AUTO','CDX']) {
    const report=JSON.parse(wasm.chemdraw_block_report_json(data,mode,c.allow_metadata_loss));
    assert.equal(report.format,'CDX');assert.deepEqual(Buffer.from(report.source_bytes),data);
    assert.equal(report.omissions.length,Number(c.allow_metadata_loss));
    assert.ok(report.omissions.every(o=>o.tag===0x0805));
    const mols=wasm.mols_from_chemdraw_block(data,mode,c.allow_metadata_loss);
    assert.equal(mols.length,report.molecule_count);
    const xml=JSON.parse(wasm.chemdraw_block_report_json(Buffer.from(c.cdxml),'CDXML',false));
    assert.deepEqual(report.coordinates,xml.coordinates);
    for (const [i,m] of mols.entries()) {
      assert.equal(m.atom_count(),report.coordinates[i].length);
      assert.match(wasm.to_mol_v3000_block(m),/V3000/);
      m.free();
    }
  }
  if(c.allow_metadata_loss)assert.throws(()=>wasm.chemdraw_block_report_json(data,'AUTO',false));
}
const data=Buffer.from(fixture.cases[0].cdx_hex,'hex');
for(let i=0;i<data.length;i++)assert.throws(()=>wasm.chemdraw_block_report_json(data.subarray(0,i),'CDX',false));
for(const [bytes,format] of [[Buffer.from('no xml'),'AUTO'],[Buffer.from('<CDXML/>'),'CDX'],[Buffer.from([255]),'CDXML'],[data,'unknown'],[Buffer.concat([data,Buffer.from([0])]),'CDX']]) {
  assert.throws(()=>wasm.chemdraw_block_report_json(bytes,format,false));
  assert.throws(()=>wasm.mols_from_chemdraw_block(bytes,format,false));
}
assert.throws(()=>wasm.chemdraw_block_report_json(new Uint8Array(1_000_001),'AUTO',false));
const bom=Buffer.from('\ufeff<CDXML><page/></CDXML>');
assert.deepEqual(Buffer.from(JSON.parse(wasm.chemdraw_block_report_json(bom,'AUTO',false)).source_bytes),bom);
console.log('ChemDraw: 56 CDX recipes, explicit/AUTO dispatch, reports/handles and malformed/resource contracts passed');
