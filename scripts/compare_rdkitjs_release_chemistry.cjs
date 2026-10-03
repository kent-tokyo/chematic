#!/usr/bin/env node
/* Compare two official RDKit.js releases on the same exposed corpus.
 * This is an oracle-version comparison, not a CheMatic parity claim.
 */

const crypto = require('node:crypto');
const fs = require('node:fs');
const path = require('node:path');
const zlib = require('node:zlib');

function options(argv) {
  const result = {};
  for (let i = 0; i < argv.length; i += 2) {
    if (!argv[i].startsWith('--') || !argv[i + 1]) throw new Error('expected --name value pairs');
    result[argv[i].slice(2)] = argv[i + 1];
  }
  for (const key of ['old-package', 'new-package', 'corpus', 'queries', 'summary', 'rows-output']) {
    if (!result[key]) throw new Error(`missing --${key}`);
  }
  return result;
}

function digest(value) {
  return crypto.createHash('sha256').update(value).digest('hex');
}

function fileRecord(filename) {
  const bytes = fs.readFileSync(filename);
  return { path: filename, bytes: bytes.length, sha256: digest(bytes) };
}

function packageRecord(directory) {
  const metadata = JSON.parse(fs.readFileSync(path.join(directory, 'package.json'), 'utf8'));
  if (metadata.name !== '@rdkit/rdkit') throw new Error(`unexpected package: ${metadata.name}`);
  return {
    name: metadata.name,
    version: metadata.version,
    package_json: fileRecord(path.join(directory, 'package.json')),
    wasm: fileRecord(path.join(directory, 'dist', 'RDKit_minimal.wasm')),
    javascript: fileRecord(path.join(directory, 'dist', 'RDKit_minimal.js')),
  };
}

async function loadRDKit(directory) {
  const initialize = require(path.resolve(directory, 'dist', 'RDKit_minimal.js'));
  return initialize({ locateFile: (name) => path.resolve(directory, 'dist', name) });
}

function atomSets(raw) {
  const matches = JSON.parse(raw);
  if (matches && !Array.isArray(matches) && Object.keys(matches).length === 0) return [];
  if (!Array.isArray(matches)) throw new Error(`unexpected substructure output: ${raw}`);
  return [...new Set(matches.map((match) => JSON.stringify([...match.atoms].sort((a, b) => a - b))))]
    .sort();
}

function inspect(module, smiles, queries) {
  const molecule = module.get_mol(smiles);
  if (!molecule) return { status: 'parse_error' };
  try {
    const canonical = molecule.get_smiles();
    const morgan = molecule.get_morgan_fp(JSON.stringify({ radius: 2, nBits: 2048 }));
    const representation = JSON.parse(molecule.get_json()).molecules[0];
    const graph = { atoms: representation.atoms, bonds: representation.bonds };
    const stereo = JSON.parse(molecule.get_stereo_tags());
    const matches = queries.map((query) => query === null ? null : atomSets(molecule.get_substruct_matches(query)));
    return { status: 'ok', canonical, morgan_sha256: digest(morgan),
      graph_sha256: digest(JSON.stringify(graph)), stereo, matches };
  } finally {
    molecule.delete();
  }
}

async function main() {
  const args = options(process.argv.slice(2));
  const corpusPath = path.resolve(args.corpus);
  const queriesPath = path.resolve(args.queries);
  const oldPath = path.resolve(args['old-package']);
  const newPath = path.resolve(args['new-package']);
  const oldMetadata = packageRecord(oldPath);
  const newMetadata = packageRecord(newPath);
  if (oldMetadata.version === newMetadata.version) throw new Error('old and new versions must differ');
  const allLines = fs.readFileSync(corpusPath, 'utf8').split(/\r?\n/).filter(Boolean);
  const limit = args.limit === undefined ? allLines.length : Number(args.limit);
  if (!Number.isSafeInteger(limit) || limit < 1 || limit > allLines.length) throw new Error('invalid --limit');
  const lines = allLines.slice(0, limit);
  const queries = JSON.parse(fs.readFileSync(queriesPath, 'utf8')).queries;
  if (!Array.isArray(queries) || queries.length === 0) throw new Error('queries must be non-empty');
  const oldRDKit = await loadRDKit(oldPath);
  const newRDKit = await loadRDKit(newPath);
  const oldQueryMols = queries.map((query) => oldRDKit.get_qmol(query));
  const newQueryMols = queries.map((query) => newRDKit.get_qmol(query));
  const oldParseFailures = oldQueryMols.flatMap((query, index) => query ? [] : [index]);
  const newParseFailures = newQueryMols.flatMap((query, index) => query ? [] : [index]);
  const counts = {
    rows: lines.length, query_count: queries.length, parse_status_differences: 0,
    canonical_spelling_differences: 0, morgan_bit_differences: 0,
    graph_json_differences: 0, cip_tag_differences: 0, cip_comparable_rows: 0,
    smarts_cells_compared: 0, smarts_cell_differences: 0,
  };
  const examples = [];
  const rowRecords = [];
  try {
    for (let index = 0; index < lines.length; index++) {
      const smiles = lines[index].trim().split(/\s+/)[0];
      const old = inspect(oldRDKit, smiles, oldQueryMols);
      const current = inspect(newRDKit, smiles, newQueryMols);
      const differences = [];
      if (old.status !== current.status) {
        counts.parse_status_differences++;
        differences.push({ operation: 'parse', old: old.status, new: current.status });
      }
      if (old.status === 'ok' && current.status === 'ok') {
        if (old.canonical !== current.canonical) {
          counts.canonical_spelling_differences++;
          differences.push({ operation: 'canonical_smiles', old: old.canonical, new: current.canonical });
        }
        if (old.morgan_sha256 !== current.morgan_sha256) {
          counts.morgan_bit_differences++;
          differences.push({ operation: 'morgan_bits', old_sha256: old.morgan_sha256, new_sha256: current.morgan_sha256 });
        }
        if (old.graph_sha256 !== current.graph_sha256) {
          counts.graph_json_differences++;
          differences.push({ operation: 'graph_json', old_sha256: old.graph_sha256, new_sha256: current.graph_sha256 });
        } else {
          counts.cip_comparable_rows++;
          if (JSON.stringify(old.stereo) !== JSON.stringify(current.stereo)) {
            counts.cip_tag_differences++;
            differences.push({ operation: 'cip_tags', old: old.stereo, new: current.stereo });
          }
        }
        for (let queryIndex = 0; queryIndex < queries.length; queryIndex++) {
          if (oldQueryMols[queryIndex] === null || newQueryMols[queryIndex] === null) continue;
          counts.smarts_cells_compared++;
          if (JSON.stringify(old.matches[queryIndex]) !== JSON.stringify(current.matches[queryIndex])) {
            counts.smarts_cell_differences++;
            differences.push({ operation: 'smarts', query_index: queryIndex, query: queries[queryIndex], old: old.matches[queryIndex], new: current.matches[queryIndex] });
          }
        }
      }
      const record = { input_index: index, smiles, old_status: old.status, new_status: current.status,
        old_canonical: old.canonical ?? null, new_canonical: current.canonical ?? null,
        old_morgan_sha256: old.morgan_sha256 ?? null, new_morgan_sha256: current.morgan_sha256 ?? null,
        old_graph_sha256: old.graph_sha256 ?? null, new_graph_sha256: current.graph_sha256 ?? null,
        old_cip_sha256: old.stereo ? digest(JSON.stringify(old.stereo)) : null,
        new_cip_sha256: current.stereo ? digest(JSON.stringify(current.stereo)) : null,
        old_smarts_sha256: old.matches ? digest(JSON.stringify(old.matches)) : null,
        new_smarts_sha256: current.matches ? digest(JSON.stringify(current.matches)) : null,
        differences };
      rowRecords.push(JSON.stringify(record));
      if (differences.length && examples.length < 30) examples.push(record);
    }
  } finally {
    oldQueryMols.forEach((query) => query?.delete());
    newQueryMols.forEach((query) => query?.delete());
  }
  const rowsBody = rowRecords.join('\n') + '\n';
  const compressedRows = zlib.gzipSync(rowsBody, { level: 9 });
  const summary = {
    schema_version: 1,
    profile: 'rdkitjs_old_new_exposed_chemistry_v1',
    generated_at_utc: new Date().toISOString(),
    node: process.version,
    old: { ...oldMetadata, runtime_version: oldRDKit.version(), query_parse_failures: oldParseFailures },
    new: { ...newMetadata, runtime_version: newRDKit.version(), query_parse_failures: newParseFailures },
    corpus: fileRecord(corpusPath), queries: fileRecord(queriesPath),
    counts, examples,
    rows_output: { path: args['rows-output'], rows: rowRecords.length,
      uncompressed_sha256: digest(rowsBody), compressed_bytes: compressedRows.length,
      compressed_sha256: digest(compressedRows) },
    caveat: 'Official RDKit.js oracle-version delta on exposed inputs; no CheMatic accuracy or speed claim.',
  };
  fs.writeFileSync(args['rows-output'], compressedRows);
  fs.writeFileSync(args.summary, JSON.stringify(summary, null, 2) + '\n');
  console.log(JSON.stringify({ old: summary.old.runtime_version, new: summary.new.runtime_version, counts }));
}

if (require.main === module) {
  main().catch((error) => { console.error(error); process.exitCode = 1; });
}

module.exports = { atomSets, options };
