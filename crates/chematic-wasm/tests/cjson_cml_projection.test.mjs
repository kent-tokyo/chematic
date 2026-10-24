import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {createRequire} from 'node:module';
import {fileURLToPath} from 'node:url';
const require=createRequire(import.meta.url);
const wasm=require(process.env.CHEMATIC_WASM_MODULE??process.argv[2]??fileURLToPath(new URL('../pkg-node/chematic_wasm.js',import.meta.url)));
const packet=JSON.parse(readFileSync(new URL('../../../validation/cjson_cml_projection/cases.json',import.meta.url)));
// JSON.stringify changes -0 to 0. Strip whitespace without parsing numbers.
function compact(s){let out='',quoted=false,escaped=false;for(const c of s){if(quoted){out+=c;if(escaped)escaped=false;else if(c==='\\')escaped=true;else if(c==='"')quoted=false;}else if(c==='"'){quoted=true;out+=c;}else if(!' \t\r\n'.includes(c))out+=c;}return out;}
for(const c of packet.valid)for(const source of [c.input,c.input.replaceAll('\n','\r\n'),compact(c.input)]){
 const output=wasm.project_cjson_to_cml(source);
 assert.equal(wasm.convert_common_format(source,'cjson','cml'),output,c.id);
 const atoms=[...output.matchAll(/<atom\s[^>]*\/>/g)];
 assert.equal(atoms.length,c.expected.numbers.length);
 if(c.expected.coordinates===null)assert.ok(!output.includes(' x3='));
 else{
  const coords=atoms.flatMap(a=>['x3','y3','z3'].map(k=>Number(a[0].match(new RegExp(k+'="([^"]+)"'))[1])));
  assert.ok(coords.every((v,i)=>Object.is(v,c.expected.coordinates.flat()[i])),c.id);
 }
}
for(const c of packet.invalid)for(const call of [()=>wasm.project_cjson_to_cml(c.input),()=>wasm.convert_common_format(c.input,'cjson','cml')])assert.throws(call,undefined,c.id);
assert.throws(()=>wasm.project_cjson_to_cml(' '.repeat(1000001)),/input bytes/);
const many=n=>JSON.stringify({atoms:{elements:{number:Array(n).fill(6)}}});
assert.equal([...wasm.project_cjson_to_cml(many(10000)).matchAll(/<atom\s/g)].length,10000);
assert.throws(()=>wasm.project_cjson_to_cml(many(10001)),/atom records/);
console.log('CJSON -> CML: 27 projections, 20 direct/common refusals and byte limit passed');
