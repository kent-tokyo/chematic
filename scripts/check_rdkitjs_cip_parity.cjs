#!/usr/bin/env node
/* Compare published CheMatic WASM CIP labels with official RDKit.js.
 * Atom/bond indices are compared only after V2000 atom symbols and bond
 * endpoints prove the two parsed molecules use the same index space.
 */

const crypto = require('node:crypto');
const fs = require('node:fs');
const path = require('node:path');
const { createRequire } = require('node:module');
const { pathToFileURL } = require('node:url');
const zlib = require('node:zlib');
const requireFromHere = createRequire(__filename);

function sha256(data) {
  return crypto.createHash('sha256').update(data).digest('hex');
}

function options(argv) {
  const result = {};
  for (let i = 0; i < argv.length; i += 2) {
    if (!argv[i]?.startsWith('--') || !argv[i + 1]) throw new Error('expected --name value pairs');
    result[argv[i].slice(2)] = argv[i + 1];
  }
  for (const name of ['rdkit-package', 'chematic-package', 'corpus', 'summary', 'rows-output']) {
    if (!result[name]) throw new Error(`missing --${name}`);
  }
  return result;
}

function pair(a, b) {
  return `${Math.min(a, b)}-${Math.max(a, b)}`;
}

function parseV2000(text) {
  const lines = text.split(/\r?\n/);
  const counts = lines.findIndex((line) => line.includes('V2000'));
  if (counts < 0) throw new Error('V2000 counts line missing');
  const atomCount = Number(lines[counts].slice(0, 3));
  const bondCount = Number(lines[counts].slice(3, 6));
  if (!Number.isSafeInteger(atomCount) || !Number.isSafeInteger(bondCount)
      || atomCount < 1 || bondCount < 0 || lines.length < counts + 1 + atomCount + bondCount) {
    throw new Error('invalid V2000 counts');
  }
  const atoms = lines.slice(counts + 1, counts + 1 + atomCount).map((line) => line.slice(31, 34).trim());
  const bonds = lines.slice(counts + 1 + atomCount, counts + 1 + atomCount + bondCount).map((line) => {
    const a = Number(line.slice(0, 3)) - 1;
    const b = Number(line.slice(3, 6)) - 1;
    const order = Number(line.slice(6, 9));
    if (!Number.isSafeInteger(a) || !Number.isSafeInteger(b)
        || a < 0 || b < 0 || a >= atomCount || b >= atomCount) throw new Error('invalid V2000 bond');
    return { a, b, order };
  });
  return { atoms, bonds };
}

function indexCorresponds(left, right) {
  return JSON.stringify(left.atoms) === JSON.stringify(right.atoms)
    && JSON.stringify(left.bonds.map(({ a, b }) => pair(a, b)).sort())
       === JSON.stringify(right.bonds.map(({ a, b }) => pair(a, b)).sort());
}

function rdkitLabels(stereo) {
  const atoms = Object.fromEntries(stereo.CIP_atoms
    .map(([index, label]) => [index, label.replace(/^\(|\)$/g, '')])
    .filter(([, label]) => ['R', 'S', 'r', 's'].includes(label)));
  const bonds = Object.fromEntries(stereo.CIP_bonds
    .map(([a, b, label]) => [pair(a, b), label.replace(/^\(|\)$/g, '')])
    .filter(([, label]) => ['E', 'Z'].includes(label)));
  return { atoms, bonds };
}

function chematicLabels(entries, graph) {
  const atoms = {};
  const bonds = {};
  for (const entry of entries) {
    const index = Number(entry.atomIdx);
    const label = String(entry.cipCode);
    if (['R', 'S', 'r', 's'].includes(label)) atoms[index] = label;
    else if (label === 'E' || label === 'Z') {
      const candidates = graph.bonds.filter(({ a, order }) => a === index && order === 2);
      if (candidates.length !== 1) throw new Error(`ambiguous E/Z bond at atom ${index}`);
      bonds[pair(candidates[0].a, candidates[0].b)] = label;
    } else throw new Error(`unknown CIP label ${label}`);
  }
  return { atoms, bonds };
}

function differingLabels(left, right) {
  const differences = [];
  for (const kind of ['atoms', 'bonds']) {
    const keys = new Set([...Object.keys(left[kind]), ...Object.keys(right[kind])]);
    for (const key of keys) {
      if (left[kind][key] !== right[kind][key]) differences.push({ kind, key, chematic: left[kind][key] ?? null, rdkit: right[kind][key] ?? null });
    }
  }
  return differences;
}

async function main() {
  const args = options(process.argv.slice(2));
  const rdkitDir = path.resolve(args['rdkit-package']);
  const chematicDir = path.resolve(args['chematic-package']);
  const corpusPath = path.resolve(args.corpus);
  const corpus = fs.readFileSync(corpusPath);
  const allSmiles = corpus.toString('utf8').split(/\r?\n/).filter(Boolean).map((line) => line.trim().split(/\s+/)[0]);
  const limit = args.limit === undefined ? allSmiles.length : Number(args.limit);
  if (!Number.isSafeInteger(limit) || limit < 1 || limit > allSmiles.length) throw new Error('invalid --limit');
  const rdkitInit = requireFromHere(path.join(rdkitDir, 'dist', 'RDKit_minimal.js'));
  const rdkit = await rdkitInit({ locateFile: (name) => path.join(rdkitDir, 'dist', name) });
  const chematic = await import(pathToFileURL(path.join(chematicDir, 'chematic_wasm.js')));
  await chematic.default({ module_or_path: new Uint8Array(fs.readFileSync(path.join(chematicDir, 'chematic_wasm_bg.wasm'))) });
  const counts = { rows: limit, graph_correspondence: 0, exact_rows: 0, typed_abstention_rows: 0, mismatch_rows: 0, unproven_rows: 0 };
  const examples = [];
  const rows = [];
  for (let index = 0; index < limit; index++) {
    const smiles = allSmiles[index];
    const record = { input_index: index, smiles, status: 'unproven', differences: [] };
    let left;
    let right;
    try {
      left = chematic.parse_smiles(smiles);
      right = rdkit.get_mol(smiles);
      if (!right) throw new Error('RDKit parse failed');
      const leftGraph = parseV2000(chematic.to_mol_block(left));
      const rightGraph = parseV2000(right.get_molblock());
      if (!indexCorresponds(leftGraph, rightGraph)) throw new Error('atom order or bond endpoints differ');
      counts.graph_correspondence++;
      const leftEntries = JSON.parse(chematic.cip_assignments_accurate_json(left));
      const unresolved = JSON.parse(chematic.cip_unresolved_json(left));
      if (!Array.isArray(leftEntries) || !Array.isArray(unresolved)) throw new Error('CheMatic CIP result not an array');
      const leftLabels = chematicLabels(leftEntries, leftGraph);
      const rightLabels = rdkitLabels(JSON.parse(right.get_stereo_tags()));
      const differences = differingLabels(leftLabels, rightLabels);
      record.differences = differences;
      record.unresolved = unresolved;
      record.chematic_labels_sha256 = sha256(JSON.stringify(leftLabels));
      record.rdkit_labels_sha256 = sha256(JSON.stringify(rightLabels));
      if (differences.length === 0 && unresolved.length === 0) {
        record.status = 'exact';
        counts.exact_rows++;
      } else if (differences.length > 0 && differences.every((difference) =>
        difference.kind === 'atoms' && unresolved.some((item) => Number(item.atomIdx) === Number(difference.key)))) {
        record.status = 'typed_abstention';
        counts.typed_abstention_rows++;
      } else if (differences.length === 0 && unresolved.length > 0) {
        record.status = 'typed_abstention';
        counts.typed_abstention_rows++;
      } else {
        record.status = 'mismatch';
        counts.mismatch_rows++;
      }
    } catch (error) {
      record.reason = String(error);
      counts.unproven_rows++;
    } finally {
      left?.free();
      right?.delete();
    }
    if (record.status !== 'exact' && examples.length < 50) examples.push(record);
    rows.push(JSON.stringify(record));
  }
  const body = rows.join('\n') + '\n';
  const compressed = zlib.gzipSync(body, { level: 9 });
  const summary = {
    schema_version: 1, profile: 'published_wasm_rdkitjs_cip_parity_v1',
    generated_at_utc: new Date().toISOString(), node: process.version,
    rdkit: { version: rdkit.version(), package: JSON.parse(fs.readFileSync(path.join(rdkitDir, 'package.json'), 'utf8')).version,
      wasm_sha256: sha256(fs.readFileSync(path.join(rdkitDir, 'dist', 'RDKit_minimal.wasm'))) },
    chematic: { version: chematic.chematic_version(), package: JSON.parse(fs.readFileSync(path.join(chematicDir, 'package.json'), 'utf8')).version,
      wasm_sha256: sha256(fs.readFileSync(path.join(chematicDir, 'chematic_wasm_bg.wasm'))) },
    corpus: { path: args.corpus, sha256: sha256(corpus) }, counts, examples,
    rows_output: { path: args['rows-output'], rows: rows.length, compressed_bytes: compressed.length,
      compressed_sha256: sha256(compressed), uncompressed_sha256: sha256(body) },
    boundary: 'Atom symbols and bond endpoints from both V2000 writers must agree before indexed CIP labels are compared. Unsupported graph exports are unproven, not mismatches or exact matches.',
  };
  fs.writeFileSync(args['rows-output'], compressed);
  fs.writeFileSync(args.summary, JSON.stringify(summary, null, 2) + '\n');
  console.log(JSON.stringify(counts));
}

if (require.main === module) {
  main().catch((error) => { console.error(error); process.exitCode = 1; });
}

module.exports = { parseV2000, indexCorresponds, differingLabels, chematicLabels, rdkitLabels };
