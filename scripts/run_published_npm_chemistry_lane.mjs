#!/usr/bin/env node
/** Full exposed 10k chemistry/310k SMARTS lane on an installed npm artifact. */

import fs from 'node:fs';
import path from 'node:path';
import crypto from 'node:crypto';
import { gunzipSync } from 'node:zlib';
import { pathToFileURL } from 'node:url';

function sha256(data) { return crypto.createHash('sha256').update(data).digest('hex'); }
function fail(message) { throw new Error(message); }
function arg(name) {
  const at = process.argv.indexOf(name);
  if (at < 0 || !process.argv[at + 1]) fail(`missing ${name}`);
  return process.argv[at + 1];
}
function jsonl(file, compressed = false) {
  const raw = fs.readFileSync(file);
  return (compressed ? gunzipSync(raw) : raw).toString('utf8').trimEnd().split('\n').map(JSON.parse);
}
function normalized(matches) {
  return [...new Set(matches.map(group => JSON.stringify([...group].sort((a, b) => a - b))))]
    .map(value => JSON.parse(value))
    .sort((left, right) => {
      for (let i = 0; i < Math.min(left.length, right.length); i++) {
        if (left[i] !== right[i]) return left[i] - right[i];
      }
      return left.length - right.length;
    });
}
function same(a, b) { return JSON.stringify(a) === JSON.stringify(b); }
function sameMap(a, b) {
  return Object.keys(a).length === Object.keys(b).length &&
    Object.entries(a).every(([key, value]) => b[key] === value);
}

const packageDir = arg('--package-dir');
const corpusPath = arg('--corpus');
const queriesPath = arg('--queries');
const referencePath = arg('--python-reference');
const oraclePath = arg('--smarts-oracle');
const rowsPath = arg('--rows-output');
const summaryPath = arg('--output');
const tarballPath = arg('--tarball');
const expectedRdkit = arg('--expected-rdkit');
const expectedVersion = arg('--expected-version');
const pkg = JSON.parse(fs.readFileSync(path.join(packageDir, 'package.json'), 'utf8'));
if (pkg.version !== expectedVersion) fail(`npm package ${pkg.version} != ${expectedVersion}`);
const wasm = await import(pathToFileURL(path.resolve(packageDir, 'chematic_wasm.js')));
wasm.initSync({ module: fs.readFileSync(path.join(packageDir, 'chematic_wasm_bg.wasm')) });

const queriesDocument = JSON.parse(fs.readFileSync(queriesPath, 'utf8'));
const queries = queriesDocument.queries;
const corpus = fs.readFileSync(corpusPath, 'utf8').trimEnd().split('\n').map(line => line.trim().split(/\s+/)[0]);
const reference = jsonl(referencePath, true);
const oracleLines = jsonl(oraclePath);
const oracleManifest = oracleLines.shift();
if (!oracleManifest._manifest || oracleManifest.rdkit_version !== expectedRdkit ||
    oracleManifest.corpus_sha256 !== sha256(fs.readFileSync(corpusPath)) ||
    oracleManifest.queries_sha256 !== sha256(fs.readFileSync(queriesPath)) ||
    corpus.length !== reference.length || corpus.length !== oracleLines.length ||
    queries.length !== oracleManifest.query_count) fail('oracle/corpus/query correspondence');

const counts = {
  input: corpus.length, completed: 0, parse_failure: 0,
  canonical_cross_binding_differences: 0, cip_exact: 0, cip_typed_refusal: 0,
  cip_wrong_confident: 0, morgan_exact: 0, morgan_typed_refusal: 0,
  morgan_wrong_confident: 0, smarts_cells: 0, smarts_differences: 0,
};
const output = fs.openSync(rowsPath, 'w');
try {
  for (let index = 0; index < corpus.length; index++) {
    const smiles = corpus[index];
    const baseline = reference[index];
    const oracle = oracleLines[index];
    if (baseline.input_index !== index || oracle.input_index !== index ||
        baseline.smiles !== smiles || oracle.smiles !== smiles) fail(`row ${index} input mismatch`);
    const row = { input_index: index, smiles, status: 'completed', differences: [] };
    let mol;
    try {
      mol = wasm.parse_smiles(smiles);
    } catch (error) {
      row.status = 'parse_failure';
      row.error = String(error);
      counts.parse_failure++;
      fs.writeSync(output, JSON.stringify(row) + '\n');
      continue;
    }
    try {
      row.canonical = mol.canonical_smiles();
      if (row.canonical !== baseline.smiles_parse_write.chematic_canonical) {
        counts.canonical_cross_binding_differences++;
        row.differences.push({ operation: 'canonical', kind: 'cross_binding' });
      }
      const assignments = JSON.parse(wasm.cip_assignments_accurate_json(mol));
      const unresolved = JSON.parse(wasm.cip_unresolved_json(mol));
      const atoms = {}, bonds = {};
      for (const { atomIdx, cipCode } of assignments) {
        if (['R', 'S', 'r', 's'].includes(cipCode)) atoms[String(atomIdx)] = cipCode;
        else if (['E', 'Z'].includes(cipCode)) {
          const candidates = [];
          for (let bondIndex = 0; bondIndex < mol.bond_count(); bondIndex++) {
            const bond = JSON.parse(wasm.get_bond_info(mol, bondIndex));
            if (bond.atomFrom === atomIdx && bond.bondOrder === 2) {
              candidates.push(`${Math.min(bond.atomFrom, bond.atomTo)}-${Math.max(bond.atomFrom, bond.atomTo)}`);
            }
          }
          bonds[candidates.length === 1 ? candidates[0] : `ambiguous-ez-atom-${atomIdx}`] = cipCode;
        }
      }
      const unresolvedMap = Object.fromEntries(unresolved.map(item => [String(item.atomIdx), item.reason]));
      const expected = baseline.cip;
      const missing = Object.keys(expected.rdkit_atoms).filter(atom => atoms[atom] !== expected.rdkit_atoms[atom]);
      const wrong = Object.keys(atoms).some(atom => atoms[atom] !== expected.rdkit_atoms[atom]) ||
        !sameMap(bonds, expected.rdkit_bonds) || missing.some(atom => !unresolvedMap[atom]);
      row.cip = { atoms, bonds, unresolved: unresolvedMap };
      if (wrong) {
        counts.cip_wrong_confident++;
        row.differences.push({ operation: 'cip', kind: 'wrong_confident_or_unaccounted' });
      } else if (missing.length) counts.cip_typed_refusal++;
      else counts.cip_exact++;

      try {
        const fingerprint = Buffer.from(wasm.rdkit_ecfp4_bitvec(mol));
        row.morgan_sha256 = sha256(fingerprint);
        if (row.morgan_sha256 === baseline.morgan.rdkit_sha256) counts.morgan_exact++;
        else {
          counts.morgan_wrong_confident++;
          row.differences.push({ operation: 'morgan', kind: 'wrong_confident' });
        }
      } catch (error) {
        row.morgan_error = String(error);
        if (/unsupported/i.test(row.morgan_error)) counts.morgan_typed_refusal++;
        else {
          counts.morgan_wrong_confident++;
          row.differences.push({ operation: 'morgan', kind: 'untyped_failure' });
        }
      }
      row.smarts_matches = [];
      for (let queryIndex = 0; queryIndex < queries.length; queryIndex++) {
        let matches;
        try { matches = normalized(JSON.parse(wasm.smarts_match_atoms(queries[queryIndex], mol))); }
        catch (error) { matches = { error: String(error) }; }
        row.smarts_matches.push(matches);
        counts.smarts_cells++;
        if (!same(matches, oracle.matches?.[queryIndex])) {
          counts.smarts_differences++;
          row.differences.push({ operation: 'smarts', query_index: queryIndex, query: queries[queryIndex],
            actual: matches, expected: oracle.matches?.[queryIndex] });
        }
      }
      counts.completed++;
    } finally { mol.free(); }
    fs.writeSync(output, JSON.stringify(row) + '\n');
  }
} finally { fs.closeSync(output); }

const summary = {
  schema: 'published-npm-chemistry-v1',
  package: { name: pkg.name, version: pkg.version, tarball_sha256: sha256(fs.readFileSync(tarballPath)),
    wasm_sha256: sha256(fs.readFileSync(path.join(packageDir, 'chematic_wasm_bg.wasm'))) },
  oracle: { rdkit_version: expectedRdkit, corpus_sha256: sha256(fs.readFileSync(corpusPath)),
    queries_sha256: sha256(fs.readFileSync(queriesPath)), smarts_all_cells_sha256: sha256(fs.readFileSync(oraclePath)),
    python_reference_sha256: sha256(fs.readFileSync(referencePath)) },
  counts,
  rows: { path: rowsPath, sha256: sha256(fs.readFileSync(rowsPath)) },
  reaction_transform: { status: 'not_exposed', api: 'run_smirks' },
  operation_matrix: { status: 'not_yet_measured', reason: '63 Python operations need an explicit npm capability mapping' },
};
fs.writeFileSync(summaryPath, JSON.stringify(summary, null, 2) + '\n');
console.log(JSON.stringify(counts));
if (counts.input !== counts.completed + counts.parse_failure || counts.smarts_cells !== counts.completed * queries.length ||
    counts.cip_wrong_confident || counts.morgan_wrong_confident) process.exitCode = 1;
