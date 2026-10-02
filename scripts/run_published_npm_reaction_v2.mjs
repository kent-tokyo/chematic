#!/usr/bin/env node
/** Run the published npm reaction-transform binding on all exposed v2 cases. */

import fs from 'node:fs';
import path from 'node:path';
import crypto from 'node:crypto';
import { pathToFileURL } from 'node:url';

function required(name) {
  const index = process.argv.indexOf(name);
  if (index < 0 || !process.argv[index + 1]) throw new Error(`missing ${name}`);
  return process.argv[index + 1];
}
const digest = bytes => crypto.createHash('sha256').update(bytes).digest('hex');
const basePath = required('--base');
const strataPath = required('--strata');
const packageDir = required('--package-dir');
const tarballPath = required('--tarball');
const outputPath = required('--output');
const baseBytes = fs.readFileSync(basePath);
const strataBytes = fs.readFileSync(strataPath);
const base = JSON.parse(baseBytes);
const strata = JSON.parse(strataBytes);
if (base.length !== 57 || strata.schema_version !== 2 ||
    strata.base_fixture_sha256 !== digest(baseBytes)) {
  throw new Error('reaction fixture manifest mismatch');
}
const cases = [...base, ...strata.cases];
if (cases.length !== 83 || new Set(cases.map(item => item.id)).size !== cases.length) {
  throw new Error('expected 83 distinct reaction cases');
}
const pkg = JSON.parse(fs.readFileSync(path.join(packageDir, 'package.json'), 'utf8'));
if (pkg.name !== '@kent-tokyo/chematic' || pkg.version !== '1.0.30') {
  throw new Error('unexpected npm package');
}
const tarballSha256 = digest(fs.readFileSync(tarballPath));
if (tarballSha256 !== 'fd427a161c11b14e6cca51e78e0c35ff05bbc52faf0539a2ddc2b7090e534201') {
  throw new Error('published npm tarball digest mismatch');
}
const wasm = await import(pathToFileURL(path.resolve(packageDir, 'chematic_wasm.js')));
wasm.initSync({ module: fs.readFileSync(path.join(packageDir, 'chematic_wasm_bg.wasm')) });
if (typeof wasm.run_reactants !== 'function') throw new Error('reaction transform not exported');
const rows = cases.map(item => {
  let payload;
  try {
    payload = wasm.run_reactants(item.smirks, item.reactants.join('|'));
  } catch (error) {
    return { id: item.id, status: 'untyped_refusal', error_type: error?.constructor?.name ?? 'Error',
      detail: String(error) };
  }
  try {
    const sets = JSON.parse(payload);
    if (!Array.isArray(sets) || !sets.every(set => Array.isArray(set) && set.every(x => typeof x === 'string'))) {
      throw new Error('reaction result shape mismatch');
    }
    return { id: item.id, status: 'products', sets, raw_product_sets: sets.length };
  } catch (error) {
    return { id: item.id, status: 'invalid_serialization', payload, detail: String(error) };
  }
});
const report = {
  schema: 'published-npm-reaction-rows/v1',
  package: { name: pkg.name, version: pkg.version, tarball_sha256: tarballSha256 },
  fixtures: { base_sha256: digest(baseBytes), strata_sha256: digest(strataBytes) },
  rows,
};
fs.writeFileSync(outputPath, JSON.stringify(report, null, 2) + '\n');
console.log(`npm reaction v2: ${rows.length} rows, ${rows.filter(row => row.status === 'untyped_refusal').length} untyped refusals, ${rows.filter(row => row.status === 'invalid_serialization').length} invalid JSON`);
