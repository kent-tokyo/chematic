#!/usr/bin/env node
/** Output-only, explicit npm slice of the published Python 63-operation matrix.
 *
 * An operation without an equivalent public API is accounted for explicitly.
 * This is not a timing benchmark. The per-row Python output archive is the
 * comparator, and both the npm tarball and corpus are SHA-256 pinned.
 */

import fs from 'node:fs';
import path from 'node:path';
import crypto from 'node:crypto';
import { gunzipSync, gzipSync } from 'node:zlib';
import { pathToFileURL } from 'node:url';

function arg(name) {
  const at = process.argv.indexOf(name);
  if (at < 0 || !process.argv[at + 1]) throw new Error(`missing ${name}`);
  return process.argv[at + 1];
}
function sha256(raw) { return crypto.createHash('sha256').update(raw).digest('hex'); }
function same(a, b) { return JSON.stringify(a) === JSON.stringify(b); }
function composition(formula) {
  if (typeof formula !== 'string') return null;
  const counts = new Map();
  const tokens = formula.matchAll(/([A-Z][a-z]?)(\d*)/g);
  let consumed = 0;
  for (const token of tokens) {
    if (token.index !== consumed) return null;
    consumed += token[0].length;
    counts.set(token[1], (counts.get(token[1]) ?? 0) + Number(token[2] || 1));
  }
  return consumed === formula.length && consumed > 0 ? [...counts].sort() : null;
}
function sameRoundedCoordinates(actual, expected) {
  return Array.isArray(actual) && Array.isArray(expected) && actual.length === expected.length &&
    actual.every((point, atom) => Array.isArray(point) && point.length === 3 &&
      Array.isArray(expected[atom]) && expected[atom].length === 3 &&
      point.every((coordinate, axis) => typeof coordinate === 'number' &&
        typeof expected[atom][axis] === 'number' && Number.isFinite(coordinate) &&
        Number.isFinite(expected[atom][axis]) &&
        Math.abs(coordinate - expected[atom][axis]) <= 0.000051));
}
function sameMeaning(op, actual, expected) {
  if (same(actual, expected)) return 'exact';
  if (!('value' in actual) || !('value' in expected)) return 'different';
  if (['qed', 'chi1v'].includes(op) &&
      typeof actual.value === 'number' && typeof expected.value === 'number' &&
      Math.abs(actual.value - expected.value) <= 1e-12 * Math.max(1, Math.abs(expected.value))) {
    return 'numeric_roundoff';
  }
  if (op === 'formula' && same(composition(actual.value), composition(expected.value)) &&
      composition(actual.value) !== null) return 'formula_spelling_only';
  if (op === 'embed_3d' && sameRoundedCoordinates(actual.value, expected.value)) return 'coordinate_roundoff';
  return 'different';
}

const packageDir = arg('--package-dir');
const tarball = arg('--tarball');
const corpusFile = arg('--corpus');
const pythonFile = arg('--python-outputs');
const outputFile = arg('--output');
const rowsFile = arg('--rows-output');
const expectedVersion = arg('--expected-version');
const expectedCorpusSha = arg('--expected-corpus-sha256');

const tarballBytes = fs.readFileSync(tarball);
const corpusBytes = fs.readFileSync(corpusFile);
if (sha256(corpusBytes) !== expectedCorpusSha) throw new Error('corpus hash mismatch');
const pkg = JSON.parse(fs.readFileSync(path.join(packageDir, 'package.json'), 'utf8'));
if (pkg.name !== '@kent-tokyo/chematic' || pkg.version !== expectedVersion) {
  throw new Error('npm artifact identity mismatch');
}
const wasm = await import(pathToFileURL(path.resolve(packageDir, 'chematic_wasm.js')));
wasm.initSync({ module: fs.readFileSync(path.join(packageDir, 'chematic_wasm_bg.wasm')) });
if (wasm.chematic_version() !== expectedVersion) throw new Error('loaded WASM version mismatch');

const smiles = corpusBytes.toString('utf8').trimEnd().split('\n').map(line => line.trim().split(/\s+/)[0]).slice(0, 5000);
const pythonBytes = fs.readFileSync(pythonFile);
const baseline = gunzipSync(pythonBytes).toString('utf8').trimEnd().split('\n').map(JSON.parse);
if (baseline.length !== 63) throw new Error('Python operation count is not 63');
const pythonMolBlocks = baseline.find(op => op.op === 'mol_block_write')?.rows;
if (!pythonMolBlocks || pythonMolBlocks.length !== 2000 ||
    pythonMolBlocks.some((row, i) => row.input_index !== i || typeof row.value !== 'string')) {
  throw new Error('pinned Python-writer MOL blocks are missing or invalid');
}

// The names are the published Python operation matrix's exact names. These
// adapters use public npm APIs only. A compatibility bundle must not use the
// native hba_count(), which has a different definition.
function productSmiles(product) {
  try { return { mol_smiles: product.canonical_smiles() }; }
  finally { product.free(); }
}
function fpBytes(bytes) { return { bytes_hex: Buffer.from(bytes).toString('hex') }; }
const popcount8 = Uint8Array.from({ length: 256 }, (_, n) => {
  let count = 0;
  for (; n; n &= n - 1) count++;
  return count;
});
let compatibleMorganTargets;
function compatibleMorganTanimotoRow(index) {
  if (!compatibleMorganTargets) {
    compatibleMorganTargets = smiles.map(smi => {
      const mol = wasm.parse_smiles(smi);
      try { return wasm.rdkit_ecfp4_bitvec(mol); }
      finally { mol.free(); }
    });
    if (compatibleMorganTargets.some(fp => fp.length !== 256)) {
      throw new Error('RDKit-compatible Morgan length is not 2048 bits');
    }
  }
  const query = compatibleMorganTargets[index];
  return compatibleMorganTargets.map(target => {
    let intersection = 0;
    let union = 0;
    for (let i = 0; i < query.length; i++) {
      intersection += popcount8[query[i] & target[i]];
      union += popcount8[query[i] | target[i]];
    }
    return union ? Math.fround(intersection / union) : 1;
  });
}
function pythonDefaultStandardize(index) {
  const result = wasm.standardize_smiles_report_json(smiles[index], false, true, true, true);
  if (result.startsWith('error:')) throw new Error(result);
  return { mol_smiles: JSON.parse(result).smiles };
}
function pythonCipLabels(mol) {
  const centers = JSON.parse(mol.assign_cip_json()).centers;
  const bonds = Array.from({ length: mol.bond_count() }, (_, index) =>
    JSON.parse(wasm.get_bond_info(mol, index)));
  return centers.map(({ atom, code }) => {
    if (code !== 'E' && code !== 'Z') return { atom_idx: atom, descriptor: code };
    const matches = bonds.flatMap((bond, index) =>
      bond.atomFrom === atom && bond.bondOrder === 2 ? [{ index, bond }] : []);
    if (matches.length !== 1) throw new Error(`ambiguous E/Z bond for atom ${atom}`);
    const { index, bond } = matches[0];
    return { atom_idx: atom, bond_atoms: [bond.atomFrom, bond.atomTo], bond_idx: index, descriptor: code };
  });
}
const smartsPatterns = ['c1ccccc1', '[OH]', 'C(=O)N', '[#7;R]', 'c1ccc2ccccc2c1',
  '[CX3](=O)[OX2H1]', '[NX3;H2,H1;!$(NC=O)]', '*~*~*~*~*~*'];
const adapters = {
  parse_smiles: (m) => ({ mol_smiles: m.canonical_smiles() }),
  'canonical_smiles(prepared)': (m) => m.canonical_smiles(),
  'parse+canonical_smiles': (m) => m.canonical_smiles(),
  mol_block_write: (m) => wasm.to_mol_block(m),
  mol_block_read: (_m, _p, index) => productSmiles(wasm.mol_from_sdf_block(pythonMolBlocks[index].value)),
  inchi: (m) => m.to_inchi(),
  inchikey: (m) => m.to_inchikey(),
  mw: (_m, p) => p.molecular_weight,
  exact_mass: (m) => m.exact_mass(),
  logp: (m) => m.logp_crippen(),
  mr: (m) => m.molar_refractivity(),
  tpsa: (m) => m.tpsa(),
  hbd: (m) => m.hbd_count(),
  hba: (_m, p) => p.hba,
  rotatable_bonds: (m) => m.rotatable_bond_count(),
  fsp3: (m) => m.fsp3(),
  ring_count: (m) => m.ring_count(),
  aromatic_ring_count: (m) => m.aromatic_ring_count(),
  qed: (m) => m.qed(),
  kappa2: (m) => m.kappa2(),
  chi1v: (m) => m.chi1v(),
  bertz_ct: (m) => m.bertz_ct(),
  formula: (m) => m.formula(),
  num_stereocenters: (m) => m.num_stereocenters(),
  'lipinski_bundle(mw,logp,hbd,hba)': (m, p) => [p.molecular_weight, m.logp_crippen(), m.hbd_count(), p.hba],
  'parse+logp': (m) => m.logp_crippen(),
  'parse+tpsa': (m) => m.tpsa(),
  'parse+ring_count': (m) => m.ring_count(),
  'parse+aromatic_ring_count': (m) => m.aromatic_ring_count(),
  'morgan_r2_2048(native ecfp4)': (m) => fpBytes(wasm.ecfp4_bitvec(m)),
  'morgan_r2_2048(rdkit-compatible)': (m) => fpBytes(wasm.rdkit_ecfp4_bitvec(m)),
  maccs: (m) => fpBytes(wasm.maccs_bitvec(m)),
  'torsion(rdkit-compatible)': (m) => fpBytes(wasm.rdkit_torsion_bitvec(m)),
  'rdkit_fp(rdkit-compatible)': (m) => fpBytes(wasm.rdkit_rdk_bitvec(m)),
  'tanimoto_1xN(per target, compatible Morgan)': (_m, _p, index) => compatibleMorganTanimotoRow(index),
  'find_matches [#6]~[#7]': (m) => JSON.parse(wasm.smarts_match_atoms('[#6]~[#7]', m)),
  murcko_scaffold: (m) => productSmiles(wasm.murcko_scaffold(m)),
  add_hydrogens: (m) => productSmiles(wasm.add_hydrogens(m)),
  remove_hydrogens: (m) => productSmiles(wasm.remove_hydrogens(m)),
  largest_fragment: (m) => productSmiles(wasm.largest_fragment(m)),
  neutralize: (m) => productSmiles(wasm.neutralize_charges(m)),
  'standardize vs Cleanup': (_m, _p, index) => pythonDefaultStandardize(index),
  'standardize vs Cleanup+Uncharge+Canonicalize': (_m, _p, index) => pythonDefaultStandardize(index),
  canonical_tautomer: (m) => productSmiles(wasm.canonical_tautomer(m)),
  brics_fragments: (m) => JSON.parse(wasm.brics_fragments_json(m)).map(s => ({ mol_smiles: s })),
  sssr_rings: (m) => {
    const membership = Array.from({ length: m.atom_count() }, () => []);
    JSON.parse(wasm.sssr_rings_json(m)).forEach((ring, ringIndex) => {
      for (const atomIndex of ring) membership[atomIndex].push(ringIndex);
    });
    return membership;
  },
  cip_labels: (m) => pythonCipLabels(m),
  svg_depiction: (m) => m.depict_svg(),
  '2d_layout': (m) => JSON.parse(wasm.depict_data_json(m)),
  embed_3d: (m) => JSON.parse(wasm.generate_3d_etkdg_coords_json(m)),
  mcs_pair: (_m, _p, index) => {
    const result = wasm.mcs_smiles_json(JSON.stringify(smiles.slice(2 * index, 2 * index + 2)));
    return result === 'null' ? null : { mol_smiles: result };
  },
};
for (const pattern of smartsPatterns) {
  adapters[`has_substructure ${pattern}`] = (m) => JSON.parse(wasm.smarts_match_atoms(pattern, m)).length > 0;
}
if (Object.keys(adapters).some(name => !baseline.some(op => op.op === name))) {
  throw new Error('adapter names drifted from Python 63-operation matrix');
}
for (const name of ['parse_smiles', 'get_rdkit_descriptors_json']) {
  if (typeof wasm[name] !== 'function') throw new Error(`npm export missing: ${name}`);
}

const selected = baseline.filter(op => Object.hasOwn(adapters, op.op));
if (selected.some(op => op.rows.length <= 0 || op.rows.length > smiles.length)) {
  throw new Error('adapted Python operation has an invalid input denominator');
}
const outputRows = selected.map(op => ({ op: op.op, rows: [] }));
for (let index = 0; index < smiles.length; index++) {
  let mol;
  try {
    mol = wasm.parse_smiles(smiles[index]);
    const profile = JSON.parse(wasm.get_rdkit_descriptors_json(mol));
    for (const [opIndex, record] of outputRows.entries()) {
      if (index >= selected[opIndex].rows.length) continue;
      try { record.rows.push({ input_index: index, value: adapters[record.op](mol, profile, index) }); }
      catch (error) { record.rows.push({ input_index: index, error: { type: 'Error', message: String(error) } }); }
    }
  } catch (error) {
    throw new Error(`row ${index} parse/profile failed: ${error}`);
  } finally {
    mol?.free();
  }
}

const comparison = outputRows.map((actual, opIndex) => {
  const expected = selected[opIndex];
  const outcomes = { exact: 0, numeric_roundoff: 0, formula_spelling_only: 0,
    coordinate_roundoff: 0, different: 0 };
  const samples = [];
  for (let i = 0; i < expected.rows.length; i++) {
    const outcome = sameMeaning(actual.op, actual.rows[i], expected.rows[i]);
    outcomes[outcome]++;
    if (outcome !== 'exact') {
      if (samples.length < 3) samples.push({ input_index: i, actual: actual.rows[i], expected: expected.rows[i] });
    }
  }
  return { op: actual.op, input_count: expected.rows.length, outcomes, samples };
});
const rawRows = Buffer.from(outputRows.map(row => JSON.stringify(row)).join('\n') + '\n');
fs.writeFileSync(rowsFile, gzipSync(rawRows, { level: 9, mtime: 0 }));
const noEquivalentPublicApi = {
  rdkit_tpsa: 'no exported N/O-only RDKit TPSA profile; MolHandle.tpsa and RDKit descriptor JSON use native TPSA',
  'atom_pair(rdkit-compatible)': 'atom_pair_bitvec is native, not the Python RDKit-compatible API',
  'pattern_fp(rdkit-compatible)': 'no matching exported bit-vector API in the published TypeScript surface',
  'embed+minimize_mmff94': 'published MMFF API generates its own coordinates and returns an energy summary, not optimized coordinates from explicit ETKDG input',
};
const notExposed = baseline.filter(op => !Object.hasOwn(adapters, op.op)).map(op => op.op);
if (!same(Object.keys(noEquivalentPublicApi).sort(), [...notExposed].sort())) {
  throw new Error('non-equivalent public API ledger drifted from adapters');
}
const summary = {
  schema: 'published-npm-python-63op-output-slice/v2',
  scope: 'output only; not a timing or RDKit parity claim',
  npm: { version: pkg.version, tarball_sha256: sha256(tarballBytes) },
  python_archive_sha256: sha256(pythonBytes),
  corpus: { sha256: sha256(corpusBytes), input_count: smiles.length },
  coverage: { total_operations: baseline.length, adapted: outputRows.length,
    no_equivalent_public_api: noEquivalentPublicApi },
  output: { rows_sha256: sha256(rawRows), compressed_sha256: sha256(fs.readFileSync(rowsFile)) },
  comparison,
};
fs.writeFileSync(outputFile, JSON.stringify(summary, null, 2) + '\n');
console.log(`npm 63-op slice: ${outputRows.length}/63 adapted, ${comparison.filter(op => op.outcomes.exact === op.input_count).length} exact operations, ${comparison.reduce((n, op) => n + op.outcomes.different, 0)} value-different rows`);
