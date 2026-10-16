import assert from 'node:assert/strict';
import { createRequire } from 'node:module';
import { fileURLToPath } from 'node:url';
const require = createRequire(import.meta.url);
const modulePath = process.argv[2]
  ?? process.env.CHEMATIC_WASM_MODULE
  ?? fileURLToPath(new URL('../pkg-node/chematic_wasm.js', import.meta.url));
const wasm = require(modulePath);
const write = (smiles, coords, charges, name='LIG') => wasm.smiles_to_pdbqt_checked(smiles, JSON.stringify(coords), JSON.stringify(charges), name);
const coords = [[0,-2.345,6.789],[1,-2.345,6.789],[2,-2.345,6.789],[3,-2.345,6.789],[4,-2.345,6.789]];
const text=write('[H]OC(Cl)Br',coords,[-.123,-.123,-.123,-.123,-.123]);
assert.equal(text,wasm.smiles_to_pdbqt('[H]OC(Cl)Br',JSON.stringify(coords),'[-0.123,-0.123,-0.123,-0.123,-0.123]','LIG'));
const lines=text.split('\n').filter(s=>s.startsWith('ATOM'));
assert.equal(lines.length,5);
for(let i=0;i<5;i++) {
  assert.equal(lines[i].length,79);
  assert.equal(Number(lines[i].slice(6,11)),i+1);
  assert.deepEqual([[30,38],[38,46],[46,54]].map(([a,b])=>Number(lines[i].slice(a,b))),coords[i]);
  assert.equal(Number(lines[i].slice(70,76)),-.123);
  assert.equal(lines[i].slice(77,79).trim(),['HD','OA','C','Cl','Br'][i]);
}
for(const bad of [2261608.729,-1000,10000,9999.9999,-999.9999]) {
  for(let axis=0;axis<3;axis++) {
    const p=[0,0,0];p[axis]=bad;
    assert.throws(()=>write('C',[p],[0]),/does not fit/);
  }
}
for(const q of [100,-10,99.9999,-9.9999]) assert.throws(()=>write('C',[[0,0,0]],[q]),/charge/);
for(const count of [0,2]) {
  assert.throws(()=>write('C',Array(count).fill([0,0,0]),[0]),/coordinates/);
  assert.throws(()=>write('C',[[0,0,0]],Array(count).fill(0)),/charges/);
}
for(const name of ['\nATOM','a\rb','\0','\t','\x7f','分子','x'.repeat(1025)]) assert.throws(()=>write('C',[[0,0,0]],[0],name));
for(const s of ['*','[B]','[Si]']) assert.throws(()=>write(s,[[0,0,0]],[0]),/AutoDock type/);
assert.throws(()=>write('',[],[]),/empty/);
assert.throws(()=>wasm.smiles_to_pdbqt_checked('C','[[NaN,0,0]]','[0]','LIG'));
assert.throws(()=>wasm.smiles_to_pdbqt_checked('C','[[0,0]]','[0]','LIG'));
assert.throws(()=>wasm.smiles_to_pdbqt_checked('C','[[0,0,0]]','[null]','LIG'));
assert.throws(()=>wasm.smiles_to_pdbqt_checked('C',' '.repeat(100001),'[0]','LIG'),/size limits/);
assert.throws(()=>wasm.smiles_to_pdbqt_checked('C','[[1e400,0,0]]','[0]','LIG'));
assert.throws(()=>write('C'.repeat(10001),[],[]),/atom count/);
for(const name of ['', 'with spaces', 'x'.repeat(1024)]) {
  const line=write('C',[[9999.9994,-999.9994,-.0001]],[99.9994],name).split('\n').find(s=>s.startsWith('ATOM'));
  assert.equal(line.length,79);
  assert.equal(line.slice(30,54),'9999.999-999.999  -0.000');
  assert.equal(line.slice(70,76),'99.999');
}
console.log('PDBQT checked Node/WASM contract passed');
