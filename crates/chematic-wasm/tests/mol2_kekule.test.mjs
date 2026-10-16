import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { createRequire } from 'node:module';
import { fileURLToPath } from 'node:url';
const require = createRequire(import.meta.url);
const modulePath = process.argv[2] ?? process.env.CHEMATIC_WASM_MODULE ?? fileURLToPath(new URL('../pkg-node/chematic_wasm.js', import.meta.url));
const wasm = require(modulePath);
const cases = JSON.parse(readFileSync(new URL('../../../validation/openbabel_issues/mol2-kekule/cases.json', import.meta.url), 'utf8')).cases;
for (const c of cases) {
  for (const newline of ['\n', '\r\n']) {
    const text = c.text.replaceAll('\r\n', '\n').replaceAll('\n', newline);
    const output = wasm.kekulize_mol2_record(text);
    const before = text.match(/[^\n]*\n|[^\n]+$/g);
    const after = output.match(/[^\n]*\n|[^\n]+$/g);
    assert.equal(before.length, after.length);
    let bonds = false, changes = 0;
    before.forEach((line, i) => {
      if (line.trim().toUpperCase().startsWith('@<TRIPOS>')) bonds = line.trim().toUpperCase() === '@<TRIPOS>BOND';
      const fields = [...line.matchAll(/\S+/g)];
      if (bonds && fields[3]?.[0] === 'ar') {
        const start = fields[3].index;
        assert.match(after[i][start], /^[12]$/);
        assert.equal(after[i], line.slice(0, start) + after[i][start] + line.slice(start + 2));
        changes++;
      } else assert.equal(line, after[i]);
    });
    assert.equal(changes, c.aromatic_bonds);
    assert.equal(wasm.kekulize_mol2_record(output), output);
    assert.equal(wasm.roundtrip_mol2_record(text), text);
  }
}
function ring(n, element = 'C') {
  const atoms = Array.from({length:n}, (_,i) => `${i+1} ${element}${i+1} 0 0 0 ${element}.ar\n`).join('');
  const bonds = Array.from({length:n}, (_,i) => `${i+1} ${i+1} ${(i+1)%n+1} ar\n`).join('');
  return `@<TRIPOS>MOLECULE\nring\n${n} ${n}\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n${atoms}@<TRIPOS>BOND\n${bonds}`;
}
const rich = readFileSync(new URL('../../../validation/openbabel_issues/mol2-kekule/rich-source.mol2', import.meta.url), 'utf8');
const richOutput = wasm.kekulize_mol2_record(rich);
assert.ok(richOutput.includes('# bond comments stay ar unchanged\r\n'));
assert.ok(richOutput.includes('FLAG_ar'));
assert.ok(richOutput.includes('café\r\n'));
assert.equal(wasm.kekulize_mol2_record(richOutput), richOutput);
for (const invalid of [ring(5), ring(4), ring(6, 'O'), ring(6, 'B'), ring(6, 'N') + '@<TRIPOS>UNITY_ATOM_ATTR\n1 1\ncharge 1\n', ring(6).replace('1 1 2 ar', '1 1 2 un'), ring(513), ring(6) + ring(6), '#' + 'x'.repeat(8*1024*1024) + '\n' + ring(6)]) {
  assert.throws(() => wasm.kekulize_mol2_record(invalid));
}
console.log(`Passed ${cases.length*2} Node/WASM token-preservation cases and 9 refusal contracts`);
