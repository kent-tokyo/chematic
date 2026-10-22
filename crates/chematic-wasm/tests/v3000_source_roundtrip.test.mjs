import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {createRequire} from 'node:module';
import {fileURLToPath} from 'node:url';
import {createHash} from 'node:crypto';
const require=createRequire(import.meta.url);
const wasm=require(process.argv[2]??process.env.CHEMATIC_WASM_MODULE??fileURLToPath(new URL('../pkg-node/chematic_wasm.js',import.meta.url)));
const packet=JSON.parse(readFileSync(new URL('../../../validation/v3000_source_roundtrip/cases.json',import.meta.url)));
for(const c of packet.cases){
 assert.equal(createHash('sha256').update(c.input).digest('hex'),c.sha256);
 for(const text of [c.input,c.input.replaceAll('\n','\r\n'),c.input.trimEnd()]){
  assert.equal(wasm.roundtrip_mol_v3000_block(text),text,c.id);
  assert.equal(wasm.convert_common_format(text,'.V3000','mol-v3000'),text,c.id);
 }
}
const source=packet.cases[0].input;
for(const change of [s=>s.replace('M  END',''),s=>s+'$$$$\n',s=>s+s,
 s=>s.replace('END CTAB','END CTAB-'),s=>s.replace('M  V30 2 O','M  V30 1 O'),
 s=>s.replace('M  V30 1 C','M  V30 0 C'),s=>s.replace('M  V30 1 1 1 2','M  V30 0 1 1 2'),
 s=>s.replace('M  V30 1 1 1 2','M  V30 1 1 1 99'),s=>s.replace('0.123456789012345','NaN'),
 s=>s.replace('M  V30 END CTAB','M  V30 END CTAB\nM  V30 BEGIN CTAB\nM  V30 END CTAB'),
 s=>s.replace('M  V30 END CTAB','M  V30 END CTAB\nM  V30 ignored')]){
 assert.throws(()=>wasm.roundtrip_mol_v3000_block(change(source)));
}
assert.throws(()=>wasm.roundtrip_mol_v3000_block(' '.repeat(1000001)),/input bytes/);
const many=n=>['count test','','',source.split('\n')[3],'M  V30 BEGIN CTAB',`M  V30 COUNTS ${n} 0 0 0 0`,'M  V30 BEGIN ATOM',
 ...Array.from({length:n},(_,i)=>`M  V30 ${i+1} C 0 0 0 0`),'M  V30 END ATOM','M  V30 BEGIN BOND','M  V30 END BOND','M  V30 END CTAB','M  END',''].join('\n');
const limit=many(10000);assert.equal(wasm.roundtrip_mol_v3000_block(limit),limit);
assert.throws(()=>wasm.roundtrip_mol_v3000_block(many(10001)),/atoms/);
console.log('V3000 source: 18 exact inputs in both APIs, invalid-source and WASM bounds passed');
