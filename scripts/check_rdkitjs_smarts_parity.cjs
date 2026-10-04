#!/usr/bin/env node
/* Complete-row SMARTS atom-set comparison of published WASM packages.
 * Atom indices are used only where the companion CIP gate has proved graph
 * correspondence for the same package pair, input row, and corpus.
 */

const crypto = require('node:crypto');
const fs = require('node:fs');
const path = require('node:path');
const { createRequire } = require('node:module');
const { pathToFileURL } = require('node:url');
const zlib = require('node:zlib');
const { atomSets } = require('./compare_rdkitjs_release_chemistry.cjs');
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
  for (const name of ['rdkit-package', 'chematic-package', 'corpus', 'queries',
    'correspondence-summary', 'correspondence-rows', 'summary', 'rows-output']) {
    if (!result[name]) throw new Error(`missing --${name}`);
  }
  return result;
}

function chematicAtomSets(raw) {
  const matches = JSON.parse(raw);
  if (!Array.isArray(matches)) throw new Error(`unexpected CheMatic match output: ${raw}`);
  return [...new Set(matches.map((atoms) => JSON.stringify([...atoms].sort((a, b) => a - b))))].sort();
}

function loadCorrespondence(summaryPath, rowsPath, corpusHash, runtimeVersion, rowCount) {
  const summary = JSON.parse(fs.readFileSync(summaryPath, 'utf8'));
  const compressed = fs.readFileSync(rowsPath);
  if (summary.corpus.sha256 !== corpusHash || summary.rdkit.version !== runtimeVersion
      || summary.counts.rows !== rowCount || summary.rows_output.compressed_sha256 !== sha256(compressed)) {
    throw new Error('correspondence packet does not match corpus and RDKit runtime');
  }
  const body = zlib.gunzipSync(compressed);
  if (summary.rows_output.uncompressed_sha256 !== sha256(body)) throw new Error('correspondence rows checksum mismatch');
  const rows = body.toString('utf8').trimEnd().split('\n').map(JSON.parse);
  if (rows.length !== rowCount || rows.some((row, index) => row.input_index !== index)) {
    throw new Error('correspondence rows incomplete or reordered');
  }
  return { summary, rows };
}

async function main() {
  const args = options(process.argv.slice(2));
  const rdkitDir = path.resolve(args['rdkit-package']);
  const chematicDir = path.resolve(args['chematic-package']);
  const corpusBytes = fs.readFileSync(args.corpus);
  const allSmiles = corpusBytes.toString('utf8').split(/\r?\n/).filter(Boolean).map((line) => line.trim().split(/\s+/)[0]);
  const limit = args.limit === undefined ? allSmiles.length : Number(args.limit);
  if (!Number.isSafeInteger(limit) || limit < 1 || limit > allSmiles.length) throw new Error('invalid --limit');
  const queriesBytes = fs.readFileSync(args.queries);
  const queries = JSON.parse(queriesBytes).queries;
  if (!Array.isArray(queries) || queries.length === 0) throw new Error('invalid queries');
  const rdkitInit = requireFromHere(path.join(rdkitDir, 'dist', 'RDKit_minimal.js'));
  const rdkit = await rdkitInit({ locateFile: (name) => path.join(rdkitDir, 'dist', name) });
  const chematic = await import(pathToFileURL(path.join(chematicDir, 'chematic_wasm.js')));
  await chematic.default({ module_or_path: new Uint8Array(fs.readFileSync(path.join(chematicDir, 'chematic_wasm_bg.wasm'))) });
  const { summary: proofSummary, rows: proofRows } = loadCorrespondence(
    args['correspondence-summary'], args['correspondence-rows'], sha256(corpusBytes), rdkit.version(), allSmiles.length,
  );
  if (proofSummary.chematic.package !== JSON.parse(fs.readFileSync(path.join(chematicDir, 'package.json'), 'utf8')).version) {
    throw new Error('correspondence CheMatic package mismatch');
  }
  const rdkitQueries = queries.map((query) => rdkit.get_qmol(query));
  const counts = { rows: limit, queries: queries.length, total_cells: limit * queries.length,
    exact_cells: 0, mismatch_cells: 0, typed_refusal_cells: 0,
    error_cells: 0, unproven_cells: 0 };
  const examples = [];
  const records = [];
  try {
    for (let index = 0; index < limit; index++) {
      const smiles = allSmiles[index];
      const record = { input_index: index, smiles, status: 'compared', differences: [] };
      if (proofRows[index].status === 'unproven') {
        record.status = 'unproven_index_correspondence';
        record.reason = proofRows[index].reason;
        counts.unproven_cells += queries.length;
      } else {
        let left;
        let right;
        try {
          left = chematic.parse_smiles(smiles);
          right = rdkit.get_mol(smiles);
          if (!right) throw new Error('RDKit parse failed despite correspondence proof');
          for (let queryIndex = 0; queryIndex < queries.length; queryIndex++) {
            const query = queries[queryIndex];
            const rdkitQuery = rdkitQueries[queryIndex];
            if (!rdkitQuery) {
              counts.error_cells++;
              record.differences.push({ query_index: queryIndex, query, status: 'rdkit_query_parse_error' });
              continue;
            }
            const expected = atomSets(right.get_substruct_matches(rdkitQuery));
            try {
              const observed = chematicAtomSets(chematic.smarts_match_atoms(query, left));
              if (JSON.stringify(observed) === JSON.stringify(expected)) counts.exact_cells++;
              else {
                counts.mismatch_cells++;
                record.differences.push({ query_index: queryIndex, query, status: 'mismatch', chematic: observed, rdkit: expected });
              }
            } catch (error) {
              const message = String(error);
              const status = /unsupported|budget|not supported/i.test(message) ? 'typed_refusal' : 'error';
              counts[status === 'typed_refusal' ? 'typed_refusal_cells' : 'error_cells']++;
              record.differences.push({ query_index: queryIndex, query, status, reason: message });
            }
          }
        } finally {
          left?.free();
          right?.delete();
        }
      }
      if ((record.differences.length || record.status !== 'compared') && examples.length < 40) examples.push(record);
      records.push(JSON.stringify(record));
    }
  } finally {
    rdkitQueries.forEach((query) => query?.delete());
  }
  const body = records.join('\n') + '\n';
  const compressed = zlib.gzipSync(body, { level: 9 });
  const summary = {
    schema_version: 1, profile: 'published_wasm_rdkitjs_smarts_atom_sets_v1',
    generated_at_utc: new Date().toISOString(), node: process.version,
    rdkit: { version: rdkit.version(), package: JSON.parse(fs.readFileSync(path.join(rdkitDir, 'package.json'), 'utf8')).version,
      wasm_sha256: sha256(fs.readFileSync(path.join(rdkitDir, 'dist', 'RDKit_minimal.wasm'))) },
    chematic: { version: chematic.chematic_version(), package: JSON.parse(fs.readFileSync(path.join(chematicDir, 'package.json'), 'utf8')).version,
      wasm_sha256: sha256(fs.readFileSync(path.join(chematicDir, 'chematic_wasm_bg.wasm'))) },
    corpus: { path: args.corpus, sha256: sha256(corpusBytes) },
    queries: { path: args.queries, sha256: sha256(queriesBytes) },
    correspondence: { summary: args['correspondence-summary'], rows: args['correspondence-rows'],
      rows_sha256: proofSummary.rows_output.uncompressed_sha256 },
    counts, examples,
    rows_output: { path: args['rows-output'], rows: records.length, compressed_bytes: compressed.length,
      compressed_sha256: sha256(compressed), uncompressed_sha256: sha256(body) },
    boundary: 'Atom-set parity on rows with CIP-packet atom order and bond endpoints proved. Parser errors, typed refusals, and unproven index correspondence remain separate cells.',
  };
  fs.writeFileSync(args['rows-output'], compressed);
  fs.writeFileSync(args.summary, JSON.stringify(summary, null, 2) + '\n');
  console.log(JSON.stringify(counts));
}

if (require.main === module) {
  main().catch((error) => { console.error(error); process.exitCode = 1; });
}

module.exports = { chematicAtomSets, loadCorrespondence };
