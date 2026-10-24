import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {createRequire} from 'node:module';
import {fileURLToPath} from 'node:url';
const require=createRequire(import.meta.url),w=require(process.argv[2]??process.env.CHEMATIC_WASM_MODULE??fileURLToPath(new URL('../pkg-node/chematic_wasm.js',import.meta.url)));
const p=JSON.parse(readFileSync(new URL('../../../validation/cml_cjson_projection/cases.json',import.meta.url)));
const variants=s=>[s,s.replaceAll('><','>\n<'),s.replaceAll('><','>\r\n<')];
for(const c of p.valid)for(const s of variants(c.input)){
 const out=w.project_cml_to_cjson(s);
 assert.equal(out,w.convert_common_format(s,'cml','cjson'),c.id);
 assert.ok(JSON.parse(out).atoms.elements.number.length,c.id);
}
for(const c of p.invalid)for(const s of variants(c.input)){
 assert.throws(()=>w.project_cml_to_cjson(s),undefined,c.id);
 assert.throws(()=>w.convert_common_format(s,'cml','cjson'),undefined,c.id);
}
console.log('CML -> CJSON checked profile passed');
