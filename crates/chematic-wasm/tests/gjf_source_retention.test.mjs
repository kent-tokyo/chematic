import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {pathToFileURL} from 'node:url';
const wasm=await import(process.env.CHEMATIC_WASM_MODULE
  ? pathToFileURL(process.env.CHEMATIC_WASM_MODULE).href
  : new URL('../pkg-node/chematic_wasm.js',import.meta.url).href);
const packet=JSON.parse(readFileSync(new URL('../../../validation/gjf-cartesian-2026-10-22/source-cases.json',import.meta.url),'utf8'));
for(const c of packet.cases){
  for(const call of [()=>wasm.roundtrip_gjf_record(c.source),()=>wasm.convert_common_format(c.source,'.GJF','com'),()=>wasm.convert_common_format(c.source,'.COM','gjf')]){
    if(c.expected==='retain')assert.deepEqual(Buffer.from(call(),'utf8'),Buffer.from(c.source,'utf8'),c.id);
    else assert.throws(call,c.id);
  }
}
assert.throws(()=>wasm.convert_common_format('C','smiles','gjf'),/same-format/);
assert.throws(()=>wasm.convert_common_format('irrelevant','gjf','smiles'),/same-format/);
const prefix='# HF/STO-3G\n\nresource control\n\n0 1\n';
const source=prefix+'H 0 0 0\n'.repeat(10_000)+'\n';
assert.equal(wasm.roundtrip_gjf_record(source),source);
assert.equal(wasm.convert_common_format(source,'gjf','com'),source);
for(const call of [()=>wasm.roundtrip_gjf_record(prefix+'H 0 0 0\n'.repeat(10_001)+'\n'),()=>wasm.convert_common_format(prefix+'H 0 0 0\n'.repeat(10_001)+'\n','com','gjf')])assert.throws(call);
const base=prefix+'C 0 0 0\n\n!';
const exact=base+'x'.repeat(1_000_000-Buffer.byteLength(base));
assert.equal(wasm.roundtrip_gjf_record(exact),exact);
assert.equal(wasm.convert_common_format(exact,'com','gjf'),exact);
assert.throws(()=>wasm.roundtrip_gjf_record(exact+'x'));
assert.throws(()=>wasm.convert_common_format(exact+'x','gjf','gjf'));
const utf8TooLarge=base+'é'.repeat(500_000);
assert.ok(utf8TooLarge.length<1_000_000);
assert.throws(()=>wasm.roundtrip_gjf_record(utf8TooLarge));
console.log('GJF exact source retention: 21 positive / 71 refused cases, aliases, BOM/EOF, byte/atom boundaries passed');
