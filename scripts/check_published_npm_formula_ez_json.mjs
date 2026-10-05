#!/usr/bin/env node
/**
 * Published-npm rerun of two v1.0.31 WASM fixes that lacked registry-package
 * evidence:
 *
 * 1. `MolHandle.formula()` uses the shared Hill-order formula: every corpus row
 *    must spell its formula exactly as the reference (published Python
 *    `Mol.formula`) does.
 * 2. `run_reactants` / `enumerate_library_2way` escape E/Z backslashes in their
 *    JSON: every response must parse, and the E/Z cases must keep a `\` or `/`.
 *
 * Usage: node check_published_npm_formula_ez_json.mjs --package-dir DIR
 *          --tarball TGZ --corpus SMI --expected-corpus-sha256 HEX
 *          --python-formulas JSON --output OUT
 */

import fs from 'node:fs';
import path from 'node:path';
import crypto from 'node:crypto';
import { pathToFileURL } from 'node:url';

function arg(name) {
  const index = process.argv.indexOf(name);
  if (index < 0 || !process.argv[index + 1]) throw new Error(`missing ${name}`);
  return process.argv[index + 1];
}
const sha256 = bytes => crypto.createHash('sha256').update(bytes).digest('hex');
const packageDir = arg('--package-dir');
const tarball = fs.readFileSync(arg('--tarball'));
const corpusBytes = fs.readFileSync(arg('--corpus'));
if (sha256(corpusBytes) !== arg('--expected-corpus-sha256')) throw new Error('corpus hash mismatch');
const reference = JSON.parse(fs.readFileSync(arg('--python-formulas')));
const pkg = JSON.parse(fs.readFileSync(path.join(packageDir, 'package.json')));
const wasm = await import(pathToFileURL(path.resolve(packageDir, 'chematic_wasm.js')));
if (typeof wasm.initSync === 'function') {
  wasm.initSync({ module: fs.readFileSync(path.join(packageDir, 'chematic_wasm_bg.wasm')) });
}

const smiles = corpusBytes.toString('utf8').trimEnd().split('\n')
  .map(line => line.trim().split(/\s+/)[0]).slice(0, 5000);
if (reference.smiles_sha256 !== sha256(Buffer.from(smiles.join('\n')))) {
  throw new Error('reference formulas were produced from different rows');
}
const formula = { rows: smiles.length, exact: 0, differ: 0, both_refused: 0, one_refused: 0, examples: [] };
smiles.forEach((text, i) => {
  let actual = null;
  try {
    const mol = wasm.parse_smiles(text);
    actual = mol.formula();
    mol.free();
  } catch { actual = null; }
  const expected = reference.formulas[i];
  if (actual === null && expected === null) formula.both_refused += 1;
  else if (actual === null || expected === null) formula.one_refused += 1;
  else if (actual === expected) formula.exact += 1;
  else {
    formula.differ += 1;
    if (formula.examples.length < 10) formula.examples.push({ row: i, smiles: text, npm: actual, python: expected });
  }
});

const ezCases = [
  ['run_reactants', '[C:1]/[C:2]=[C:3]/[C:4]>>[C:1]/[C:2]=[C:3]/[C:4]', ['C/C=C/C']],
  ['run_reactants', '[C:1]/[C:2]=[C:3]\\[C:4]>>[C:1]/[C:2]=[C:3]\\[C:4]', ['C/C=C\\C']],
  ['run_reactants', '[OH:1]>>[O:1]C', ['F/C=C\\CO']],
  ['run_reactants', '[Cl:1]>>[Br:1]', ['Cl/C=C/C=C\\Cl']],
  ['enumerate_library_2way', '[C:1][Cl].[N:2]>>[C:1][N:2]', ['F/C=C\\CCl', 'NC/C=C/F']],
];
const ez = [];
for (const [api, template, inputs] of ezCases) {
  let text;
  try {
    text = api === 'run_reactants'
      ? wasm.run_reactants(template, inputs.join('|'))
      : wasm.enumerate_library_2way(template, inputs[0], inputs[1]);
  } catch (error) {
    ez.push({ api, template, inputs, status: 'error', detail: String(error) });
    continue;
  }
  let parsed;
  try { parsed = JSON.parse(text); } catch (error) {
    ez.push({ api, template, inputs, status: 'invalid_json', detail: String(error) });
    continue;
  }
  const flat = parsed.flat(2).filter(item => typeof item === 'string');
  ez.push({ api, template, inputs, status: flat.length === 0 ? 'no_products'
    : flat.every(s => s.includes('/') || s.includes('\\')) ? 'ok' : 'stereo_lost', products: flat });
}

const report = {
  schema: 'published-npm-formula-ez-json/v1',
  artifact: { package: pkg.name, version: pkg.version, tarball_sha256: sha256(tarball),
    tarball_sha512: crypto.createHash('sha512').update(tarball).digest('base64'), node: process.version },
  corpus_sha256: sha256(corpusBytes),
  reference: reference.artifact,
  formula,
  ez_json: { cases: ez.length, ok: ez.filter(c => c.status === 'ok').length, rows: ez },
};
fs.writeFileSync(arg('--output'), JSON.stringify(report, null, 2) + '\n');
console.log(`formula ${formula.exact}/${formula.rows} exact, ${formula.differ} differ, ` +
  `${formula.one_refused} one-sided refusals; E/Z JSON ${report.ez_json.ok}/${ez.length} ok`);
process.exitCode = formula.differ || formula.one_refused || report.ez_json.ok !== ez.length ? 1 : 0;
