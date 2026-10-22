import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {pathToFileURL} from 'node:url';
const wasm=await import(process.env.CHEMATIC_WASM_MODULE?pathToFileURL(process.env.CHEMATIC_WASM_MODULE).href:new URL('../pkg-node/chematic_wasm.js',import.meta.url).href);
const packet=JSON.parse(readFileSync(new URL('../../../validation/pqr-checked-2026-10-22/cases.json',import.meta.url),'utf8'));
const bits=value=>{const b=Buffer.alloc(8);b.writeDoubleBE(value);return b.toString('hex');};
for(const c of packet.cases){
  if(c.expected==='refuse'){assert.throws(()=>wasm.write_pqr_json_checked(c.records_json),c.id);continue;}
  // Pass raw JSON. Re-stringifying parsed JS values would lose -0 and wide IDs.
  const text=wasm.write_pqr_json_checked(c.records_json);
  const lines=text.trimEnd().split('\n');assert.equal(lines.length,c.fields.length);
  for(const [i,line] of lines.entries()){
    const tokens=line.split(/\s+/);const expected=c.fields[i];
    assert.ok([10,11].includes(tokens.length));
    assert.equal(tokens[0],expected.group_pdb);assert.equal(BigInt(tokens[1]).toString(),expected.serial);
    assert.equal(tokens[2],expected.atom_name);assert.equal(tokens[3],expected.res_name);
    const hasChain=tokens.length===11;
    assert.equal(hasChain?tokens[4]:null,expected.chain_id);
    const residue=tokens[hasChain?5:4].match(/^([+-]?\d+)([A-Za-z]?)$/);assert.ok(residue);
    assert.equal(BigInt(residue[1]).toString(),expected.res_seq);assert.equal(residue[2]||null,expected.icode);
    const numeric=tokens.slice(hasChain?6:5);
    for(const [axis,key] of ['x','y','z','charge','radius'].entries())assert.equal(bits(Number(numeric[axis])),expected.bits[key],`${c.id} ${key}`);
  }
}
for(const raw of ['null','{}','[null]','[]','[true]','[{}]'])assert.throws(()=>wasm.write_pqr_json_checked(raw));
assert.throws(()=>wasm.write_pqr_json_checked(' '.repeat(1_000_001)));
assert.throws(()=>wasm.write_pqr_json_checked('é'.repeat(500_001)));
// Include the JSON transport itself at the exact byte cap.
const base=packet.cases[0].records_json;
const exact=base+' '.repeat(1_000_000-Buffer.byteLength(base));
assert.equal(wasm.write_pqr_json_checked(exact),wasm.write_pqr_json_checked(base));
assert.throws(()=>wasm.write_pqr_json_checked(exact+' '));
console.log('PQR checked output: 26 bit-exact atom-field controls, 43 refusals, strict JSON and UTF-8 byte limits passed');
