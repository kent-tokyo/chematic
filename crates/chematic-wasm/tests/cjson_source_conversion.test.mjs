import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { createRequire } from 'node:module';
import { fileURLToPath } from 'node:url';
import { createHash } from 'node:crypto';
const require = createRequire(import.meta.url);
const wasm = require(process.argv[2] ?? process.env.CHEMATIC_WASM_MODULE
  ?? fileURLToPath(new URL('../pkg-node/chematic_wasm.js', import.meta.url)));
const packet = JSON.parse(readFileSync(new URL('../../../validation/cjson_source_conversion/cases.json', import.meta.url)));
const compact = text => {
  let quoted=false, escaped=false, result='';
  for (const ch of text) {
    if (quoted) {
      result+=ch;
      if (escaped) escaped=false;
      else if (ch==='\\') escaped=true;
      else if (ch==='"') quoted=false;
    } else if (ch==='"') { quoted=true; result+=ch; }
    else if (!/[ \n\r\t]/.test(ch)) result+=ch;
  }
  return result;
};
for (const c of packet.cases) {
  assert.equal(createHash('sha256').update(c.input).digest('hex'), c.sha256);
  // Keep large integers and decimal spellings intact; JS parse/stringify is
  // intentionally excluded from the expected source-retention oracle.
  for (const text of [c.input, c.input.replaceAll('\n', '\r\n'), compact(c.input)]) {
    assert.equal(wasm.convert_common_format(text, '.CJSON', 'cjson'), text, c.id);
  }
}
const base = JSON.parse(packet.cases[0].input);
for (const edit of [v=>v.chemicalJson=0, v=>v.atoms.formalCharges=[0],
                    v=>v.atoms.coords['3d']=[0], v=>v.partialCharges={method:[0]},
                    v=>v.bonds.order=[2.2], v=>v.bonds.connections.index=[0,0],
                    v=>v.bonds.connections.index=[0,2]]) {
  const v=structuredClone(base); edit(v);
  assert.throws(()=>wasm.convert_common_format(JSON.stringify(v),'cjson','cjson'));
}
for (const suffix of [',"chemicalJson":1', ',"opaque":{"id":1,"id":2}',
                      ',"opaque":{"id":1,"\\u0069d":2}']) {
  const text=packet.cases[0].input.trimEnd().slice(0,-1)+suffix+'}';
  assert.throws(()=>wasm.convert_common_format(text,'cjson','cjson'), /duplicate object key/);
}
assert.throws(()=>wasm.convert_common_format(' '.repeat(1000001),'cjson','cjson'), /maximum input size/);
const atoms=n=>JSON.stringify({chemicalJson:1,atoms:{elements:{number:Array(n).fill(6)}}});
const atLimit=atoms(10000);
assert.equal(wasm.convert_common_format(atLimit,'cjson','cjson'),atLimit);
assert.throws(()=>wasm.convert_common_format(atoms(10001),'cjson','cjson'), /atom/);
console.log('CJSON source conversion: 18 exact LF/CRLF/compact inputs, local profile, duplicates and WASM limits passed');
