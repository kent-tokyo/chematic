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
 const cml=wasm.convert_common_format(source,'.CJSON','cml');
 // CML declares total H counts; the return leg must no longer discard them.
 assert.throws(()=>wasm.convert_common_format(cml,'cml','cjson'),/cannot preserve explicit hydrogen count/);
 // Observe the checked input graph through MolJSON, without a lossy return leg.
 const projected=JSON.parse(wasm.convert_common_format(source,'cjson','moljson'));
 const symbols={1:'H',6:'C',7:'N',8:'O',9:'F',11:'Na',17:'Cl'};
 assert.deepEqual(projected.atoms.map(a=>a.element),c.expected.numbers.map(n=>symbols[n]),c.id);
 const ids=new Map(projected.atoms.map((a,i)=>[a.id,i]));
 assert.deepEqual(projected.atoms.map(a=>a.charge),c.expected.charges,c.id);
 assert.deepEqual(projected.atoms.map(a=>a.isotope??0),c.expected.isotopes,c.id);
 assert.deepEqual(projected.bonds.flatMap(b=>[ids.get(b.source_id),ids.get(b.target_id)]),c.expected.indices,c.id);
 assert.deepEqual(projected.bonds.map(b=>b.order),c.expected.orders,c.id);
 assert.equal(wasm.convert_common_format(source,'.CJSON','cjson'),source,c.id);
}
for(const c of packet.invalid)for(const source of variants(c))for(const target of ['smiles','mol','mol_v3000','mol2','cml','cjson','moljson'])assert.throws(()=>wasm.convert_common_format(source,'cjson',target),undefined,c.id);
const many=n=>JSON.stringify({atoms:{elements:{number:Array(n).fill(6)}}});
assert.equal(JSON.parse(wasm.convert_common_format(many(10000),'cjson','cjson')).atoms.elements.number.length,10000);
assert.throws(()=>wasm.convert_common_format(many(10001),'cjson','cjson'),/atom records/);
assert.throws(()=>wasm.convert_common_format(' '.repeat(1000001),'cjson','cjson'),/input size/);
console.log('CJSON checked: 45 input graphs/source retention and CML return-leg refusals, 105 invalid inputs across seven outputs, WASM byte/atom caps passed');
