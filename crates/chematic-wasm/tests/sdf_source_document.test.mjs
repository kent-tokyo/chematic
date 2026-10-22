import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {createRequire} from 'node:module';
import {fileURLToPath} from 'node:url';
const require=createRequire(import.meta.url);
const wasm=require(process.argv[2]??process.env.CHEMATIC_WASM_MODULE??fileURLToPath(new URL('../pkg-node/chematic_wasm.js',import.meta.url)));
const packet=JSON.parse(readFileSync(new URL('../../../validation/sdf_source_document/cases.json',import.meta.url)));
for(const c of packet.cases){
 for(const source of [c.input,c.input.replaceAll('\n','\r\n'),c.input.trimEnd()]){
  assert.equal(wasm.roundtrip_sdf_source(source),source,c.id);
  assert.equal(wasm.convert_common_format(source,'.SDF','sdf'),source,c.id);
  if(c.records>1){
   for(const input of ['mol','sdf'])for(const output of ['smiles','cjson','mol','v3000']){
    assert.throws(()=>wasm.convert_common_format(source,input,output),/records.*limit 1/,c.id);
   }
  }else assert.equal(wasm.convert_common_format(source,'sdf','smiles'),'O');
 }
}
const source=packet.cases[0].input;
for(const bad of ['', ' \n', '$$$$\n', '$$$$\n'+source, source+'junk\n', source.replaceAll('M  END','M  END suffix'), source.replace('    2.1250 O','       NaN O')]){
 assert.throws(()=>wasm.roundtrip_sdf_source(bad));
 assert.throws(()=>wasm.convert_common_format(bad,'sdf','sdf'));
}
const framed=wasm.convert_common_format('CO','smiles','sdf');
assert.ok(framed.endsWith('$$$$\n'));
assert.ok(wasm.convert_common_format(framed,'sdf','smiles'));
assert.throws(()=>wasm.roundtrip_sdf_source(' '.repeat(1000001)),/input bytes/);
const many=n=>['count test','','',packet.cases[2].input.split('\n')[3],'M  V30 BEGIN CTAB',`M  V30 COUNTS ${n} 0 0 0 0`,'M  V30 BEGIN ATOM',...Array.from({length:n},(_,i)=>`M  V30 ${i+1} C 0 0 0 0`),'M  V30 END ATOM','M  V30 BEGIN BOND','M  V30 END BOND','M  V30 END CTAB','M  END','$$$$',''].join('\n');
const exact=many(5000)+many(5000);assert.equal(wasm.roundtrip_sdf_source(exact),exact);
assert.throws(()=>wasm.roundtrip_sdf_source(many(5000)+many(5001)),/total atoms/);
assert.throws(()=>wasm.roundtrip_sdf_source(many(10001)),/atoms per record/);
const zero='zero\n\n\n  0  0  0  0  0  0  0  0  0  0  0999 V2000\nM  END\n>  <NOTE>\nmetadata only\n\n$$$$\n';
assert.equal(wasm.roundtrip_sdf_source(zero),zero);
assert.throws(()=>wasm.convert_common_format(zero+zero,'sdf','smiles'),/records.*limit 1/);
const utf8Header=packet.cases[2].input.replace('original high precision','original\u2028high precision');
assert.equal(wasm.convert_common_format(utf8Header,'sdf','smiles'),'O');
console.log('SDF source: 24 whole documents exact; multi-record single-molecule refusal, invalid input, framing, zero atoms, UTF-8 headers and total/per-record atom limits passed');
