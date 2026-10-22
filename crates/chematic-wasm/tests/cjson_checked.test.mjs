import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {createRequire} from 'node:module';
import {fileURLToPath} from 'node:url';
const require=createRequire(import.meta.url);
const wasm=require(process.argv[2]??process.env.CHEMATIC_WASM_MODULE??fileURLToPath(new URL('../pkg-node/chematic_wasm.js',import.meta.url)));
const packet=JSON.parse(readFileSync(new URL('../../../validation/cjson_checked/cases.json',import.meta.url)));
function compact(s){let result='',quoted=false,escaped=false;for(const c of s){if(quoted){result+=c;if(escaped)escaped=false;else if(c==='\\')escaped=true;else if(c==='"')quoted=false;}else if(c==='"'){quoted=true;result+=c;}else if(!' \t\r\n'.includes(c))result+=c;}return result;}
const variants=c=>[c.input,c.input.replaceAll('\n','\r\n'),compact(c.input)];
for(const c of [...packet.valid,...packet.compatibility])for(const source of variants(c)){
 const data=JSON.parse(wasm.convert_common_format(source,'.CJSON','cjson')),n=c.expected.numbers.length;
 assert.deepEqual(data.atoms.elements.number,c.expected.numbers,c.id);
 assert.deepEqual(data.atoms.formalCharges??Array(n).fill(0),c.expected.charges,c.id);
 assert.deepEqual(data.atoms.isotopes??Array(n).fill(0),c.expected.isotopes,c.id);
 assert.deepEqual(data.bonds.connections.index,c.expected.indices,c.id);
 assert.deepEqual(data.bonds.order,c.expected.orders,c.id);
}
for(const c of packet.invalid)for(const source of variants(c))for(const target of ['smiles','mol','mol_v3000','mol2','cml','cjson','moljson'])assert.throws(()=>wasm.convert_common_format(source,'cjson',target),undefined,c.id);
const many=n=>JSON.stringify({atoms:{elements:{number:Array(n).fill(6)}}});
assert.equal(JSON.parse(wasm.convert_common_format(many(10000),'cjson','cjson')).atoms.elements.number.length,10000);
assert.throws(()=>wasm.convert_common_format(many(10001),'cjson','cjson'),/atom records/);
assert.throws(()=>wasm.convert_common_format(' '.repeat(1000001),'cjson','cjson'),/input size/);
console.log('CJSON checked: 18 common + 27 compatibility graphs, 105 invalid inputs across seven outputs, WASM byte/atom caps passed');
