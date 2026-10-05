#!/usr/bin/env node
/**
 * Emit checked WASM reaction rows (`run_reactants_checked(..., true)`) for the
 * 83 fixtures in the schema read by
 * `reaction_python_checked_provenance_gate.py --candidates`.
 *
 * Usage: node run_published_npm_checked_candidates.mjs --base B --strata S
 *          --package-dir DIR --tarball TGZ [--published-from URL] --output OUT
 */

import fs from 'node:fs';
import path from 'node:path';
import crypto from 'node:crypto';
import { pathToFileURL } from 'node:url';

function arg(name, required = true) {
  const index = process.argv.indexOf(name);
  if (index < 0 || !process.argv[index + 1]) {
    if (required) throw new Error(`missing ${name}`);
    return null;
  }
  return process.argv[index + 1];
}
const sha256 = bytes => crypto.createHash('sha256').update(bytes).digest('hex');
const baseBytes = fs.readFileSync(arg('--base'));
const strataBytes = fs.readFileSync(arg('--strata'));
const base = JSON.parse(baseBytes);
const strata = JSON.parse(strataBytes);
if (base.length !== 57 || strata.schema_version !== 2 ||
    strata.base_fixture_sha256 !== sha256(baseBytes)) {
  throw new Error('reaction fixture manifest mismatch');
}
const cases = [...base, ...strata.cases];
if (cases.length !== 83 || new Set(cases.map(item => item.id)).size !== cases.length) {
  throw new Error('expected 83 distinct reaction cases');
}
const packageDir = arg('--package-dir');
const tarball = fs.readFileSync(arg('--tarball'));
const pkg = JSON.parse(fs.readFileSync(path.join(packageDir, 'package.json')));
const wasmBytes = fs.readFileSync(path.join(packageDir, 'chematic_wasm_bg.wasm'));
const wasm = await import(pathToFileURL(path.resolve(packageDir, 'chematic_wasm.js')));
if (typeof wasm.initSync === 'function') wasm.initSync({ module: wasmBytes });
if (typeof wasm.run_reactants_checked !== 'function') {
  throw new Error('run_reactants_checked is not exported');
}
const DIAGNOSTICS = ['accepted_matches', 'applied_products', 'valence_rejected_matches',
  'truncated_matches'];

const rows = {};
for (const item of cases) {
  const out = JSON.parse(wasm.run_reactants_checked(item.smirks, item.reactants.join('|'), true));
  const diagnostics = Object.fromEntries(DIAGNOSTICS.map(key => [key, out[key]]));
  if (out.status !== 'products' && out.status !== 'no_match') {
    rows[item.id] = { status: out.status, reason: out.reason, detail: out.detail, diagnostics };
    continue;
  }
  const sets = out.products.map((set, i) => set.map((smiles, j) => ({
    smiles,
    // WASM spells a source as {reactant, atom}; Python as [reactant, atom].
    atom_sources: out.product_atom_sources[i][j].map(
      source => source === null ? null : [source.reactant, source.atom]),
    template_map_numbers: out.product_template_maps[i][j],
  })));
  rows[item.id] = { status: out.status, sets, diagnostics };
}
const publishedFrom = arg('--published-from', false);
const report = {
  schema: 'python-checked-reaction-candidates/v1',
  fixtures: { base_sha256: sha256(baseBytes), strata_sha256: sha256(strataBytes) },
  artifact: {
    kind: publishedFrom ? 'published_npm_tarball' : 'npm_tarball',
    published: Boolean(publishedFrom),
    published_from: publishedFrom,
    chematic_version: pkg.version,
    package: pkg.name,
    node: process.version,
    tarball_sha256: sha256(tarball),
    tarball_sha512: crypto.createHash('sha512').update(tarball).digest('base64'),
    wasm_sha256: sha256(wasmBytes),
  },
  rows,
};
const sorted = value => Array.isArray(value) ? value.map(sorted)
  : value && typeof value === 'object'
    ? Object.fromEntries(Object.keys(value).sort().map(k => [k, sorted(value[k])])) : value;
fs.writeFileSync(arg('--output'), JSON.stringify(sorted(report), null, 2) + '\n');
console.log(`${arg('--output')}: ${Object.keys(rows).length} rows`);
