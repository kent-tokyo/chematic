import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {createRequire} from 'node:module';
import {fileURLToPath} from 'node:url';
const require=createRequire(import.meta.url);
const wasm=require(process.env.CHEMATIC_WASM_MODULE??process.argv[2]??fileURLToPath(new URL('../pkg-node/chematic_wasm.js',import.meta.url)));
const packet=JSON.parse(readFileSync(new URL('../../../validation/cml_source_preservation/cases.json',import.meta.url))).cases;
for(const c of packet){
 assert.equal(wasm.roundtrip_cml_source(c.source),c.source,c.id);
 assert.equal(wasm.convert_common_format(c.source,'.CML','cml'),c.source,c.id);
}
const ns='http://www.xml-cml.org/schema';
for(const text of [`<molecule xmlns="${ns}"/>`, `<molecule xmlns="${ns}"><bond atomRefs2="missing other" order="opaque"/></molecule>`, `<molecule xmlns="${ns}" xmlns:v="urn:vendor"><v:molecule/></molecule>`]){
 assert.equal(wasm.convert_common_format(text,'cml','cml'),text);
}
for(const text of ['<molecule/>',`<cml xmlns="${ns}"><molecule/></cml>`,`<molecule xmlns="${ns}"><molecule/></molecule>`,`<!DOCTYPE molecule><molecule xmlns="${ns}"/>`,`<?xml version="1.1"?><molecule xmlns="${ns}"/>`,`<molecule xmlns="${ns}"><atom id="a" id="b"/></molecule>`,`<molecule xmlns="${ns}">&#0;</molecule>`,`<molecule xmlns="${ns}"><atom/></molecule><bad/>`]){
 for(const call of [()=>wasm.roundtrip_cml_source(text),()=>wasm.convert_common_format(text,'cml','cml')])assert.throws(call);
}
const many=n=>`<molecule xmlns="${ns}"><atomArray>`+Array.from({length:n},(_,i)=>`<atom id="a${i}" elementType="C"/>`).join('')+'</atomArray></molecule>';
assert.equal(wasm.roundtrip_cml_source(many(10000)),many(10000));
assert.throws(()=>wasm.roundtrip_cml_source(many(10001)),/atom\/bond elements/);
assert.throws(()=>wasm.roundtrip_cml_source(' '.repeat(1000001)),/input bytes/);
assert.throws(()=>wasm.convert_common_format(' '.repeat(1000001),'cml','cml'),/input size/);
console.log('CML common source retention: exact source, opaque scope, malformed XML and WASM atom/byte boundaries passed');
