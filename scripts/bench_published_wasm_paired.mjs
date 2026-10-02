#!/usr/bin/env node
/**
 * Published npm/WASM parse and compatible-Morgan gate.
 *
 * Each child loads exactly one pinned artifact. Parent preflights every input
 * for atom-count and Morgan bit identity, then runs 20 ABBA/BAAB blocks in
 * fresh Node processes. RSS is a whole-process observation, not an allocation
 * attributed to one chemistry library. This is Node evidence, not a browser
 * performance claim.
 */

import { createHash } from 'node:crypto';
import { spawnSync } from 'node:child_process';
import { readFileSync, writeFileSync } from 'node:fs';
import { dirname, join, resolve } from 'node:path';
import { fileURLToPath, pathToFileURL } from 'node:url';
import { performance } from 'node:perf_hooks';

const ROOT = resolve(dirname(fileURLToPath(import.meta.url)), '..');
const SELF = fileURLToPath(import.meta.url);
const EXPECTED_CORPUS_SHA256 = 'd6f2ba3f128296f935007f0b0813aa97b6ebcc2457e014ddca2213ddd655276c';
const OPTIONS = {
  v29: {
    kind: 'chematic', version: '1.0.29',
    tarball: 'd0af78a8d6b711a13b63985d4b36079e245441e99f1712532555dc69f11569fa',
    js: '4ca57d6e9b9c8cc703d854d0265550f5ef355d52c2bf56a1a7b8bfcc09f21eb4',
    wasm: '78d3a205d5b4a8914ff3fc9dab3832a5b668e2b886ecda5dc5c566a48b2ba919',
  },
  v30: {
    kind: 'chematic', version: '1.0.30',
    tarball: 'fd427a161c11b14e6cca51e78e0c35ff05bbc52faf0539a2ddc2b7090e534201',
    js: '4ca57d6e9b9c8cc703d854d0265550f5ef355d52c2bf56a1a7b8bfcc09f21eb4',
    wasm: 'e0a505c33a282724114e279201dae199b2e5dbe68c35c0ebaa474693a6eb4ae6',
  },
  rdkit: {
    kind: 'rdkit', version: '2026.03.6',
    tarball: '3b86b72775394ae997fb96ca854c6e7c5840927d3e899c1ef117d42bca573c87',
    js: '25ea1529b8eeb106ac1f3acb7c5d338a90533248c6223614574fefd14fee772f',
    wasm: 'f5981af4cf5bdcef861ffdf7410330bf9887f1195e39fbf9901b3a2f134ace2e',
  },
};
const LANES = ['parse_only', 'parse_morgan', 'prepared_first_use', 'prepared_reused'];
const BLOCKS = 20;
const MORGAN_OPTIONS = '{"radius":2,"nBits":2048}';

function arg(name) {
  const index = process.argv.indexOf(name);
  if (index < 0 || !process.argv[index + 1]) throw new Error(`missing ${name}`);
  return process.argv[index + 1];
}
function sha(bytes) { return createHash('sha256').update(bytes).digest('hex'); }
function fileSha(path) { return sha(readFileSync(path)); }
function check(condition, message) { if (!condition) throw new Error(message); }
function q(values, fraction) {
  const sorted = [...values].sort((a, b) => a - b);
  const index = (sorted.length - 1) * fraction;
  const low = Math.floor(index);
  return sorted[low] + (sorted[Math.min(low + 1, sorted.length - 1)] - sorted[low]) * (index - low);
}
function optionsFor(name) {
  check(Object.hasOwn(OPTIONS, name), `unknown package ${name}`);
  return OPTIONS[name];
}
function artifactPaths(name) {
  const dir = resolve(arg(`--${name}-dir`));
  const info = optionsFor(name);
  const js = join(dir, info.kind === 'rdkit' ? 'dist/RDKit_minimal.js' : 'chematic_wasm.js');
  const wasm = join(dir, info.kind === 'rdkit' ? 'dist/RDKit_minimal.wasm' : 'chematic_wasm_bg.wasm');
  const pkg = JSON.parse(readFileSync(join(dir, 'package.json'), 'utf8'));
  check(pkg.version === info.version, `${name} package version mismatch`);
  check(fileSha(js) === info.js && fileSha(wasm) === info.wasm, `${name} JS/WASM hash mismatch`);
  return { dir, js, wasm, info };
}
function rows() {
  const bytes = readFileSync(resolve(arg('--corpus')));
  check(sha(bytes) === EXPECTED_CORPUS_SHA256, 'pinned corpus SHA-256 mismatch');
  const count = Number(arg('--rows'));
  check(Number.isInteger(count) && count > 0, 'rows must be positive');
  const values = bytes.toString('utf8').split(/\r?\n/).map(x => x.trim()).filter(Boolean).slice(0, count);
  check(values.length === count && values.every(x => !/\s/.test(x)), 'corpus row contract failed');
  return { values, hash: sha(bytes), count };
}
function packedBits(bytes) {
  check(bytes.length === 256, 'CheMatic Morgan must be 2048 bits');
  let bits = '';
  for (const byte of bytes) for (let bit = 0; bit < 8; bit++) bits += (byte >> bit) & 1;
  return bits;
}
async function load(name) {
  const { js, wasm, info } = artifactPaths(name);
  const start = performance.now();
  const module = await import(pathToFileURL(js));
  if (info.kind === 'chematic') {
    module.initSync({ module: readFileSync(wasm) });
    check(module.chematic_version() === info.version, 'loaded CheMatic version mismatch');
    return {
      init_ms: performance.now() - start,
      parse: smi => {
        const mol = module.parse_smiles(smi);
        check(mol, 'CheMatic parse returned null');
        return mol;
      },
      count: mol => mol.atom_count(),
      fp: mol => packedBits(module.rdkit_ecfp4_bitvec(mol)),
      free: mol => mol.free(),
    };
  }
  const initialize = module.default ?? module.initRDKitModule;
  check(typeof initialize === 'function', 'RDKit initializer missing');
  const rdkit = await initialize({ locateFile: name => name.endsWith('.wasm') ? wasm : name });
  check(rdkit.version() === info.version, 'loaded RDKit version mismatch');
  return {
    init_ms: performance.now() - start,
    parse: smi => {
      const mol = rdkit.get_mol(smi);
      check(mol, 'RDKit parse returned null');
      return mol;
    },
    count: mol => mol.get_num_atoms(),
    fp: mol => {
      const bits = mol.get_morgan_fp(MORGAN_OPTIONS);
      check(typeof bits === 'string' && /^[01]{2048}$/.test(bits), 'RDKit Morgan format changed');
      return bits;
    },
    free: mol => mol.delete(),
  };
}
function withMol(api, text, fn) {
  const mol = api.parse(text);
  try { return fn(mol); } finally { api.free(mol); }
}
function preflight(api, values) {
  const bits = createHash('sha256');
  const rows = values.map(text => withMol(api, text, mol => {
    const fingerprint = api.fp(mol);
    bits.update(fingerprint); bits.update('\n');
    return [api.count(mol), sha(fingerprint)];
  }));
  return { rows, morgan_rows_sha256: bits.digest('hex') };
}
function warm(api, values, lane) {
  for (const text of values.slice(0, Math.min(20, values.length))) {
    withMol(api, text, mol => {
      if (lane !== 'parse_only') api.fp(mol);
    });
  }
}
function timed(api, values, lane) {
  const prepared = lane.startsWith('prepared') ? values.map(api.parse) : [];
  const hash = createHash('sha256');
  try {
    if (lane === 'prepared_reused') for (const mol of prepared) api.fp(mol);
    const start = process.hrtime.bigint();
    for (let i = 0; i < values.length; i++) {
      if (lane === 'parse_only') {
        withMol(api, values[i], () => {});
      } else if (lane === 'parse_morgan') {
        withMol(api, values[i], mol => {
          hash.update(api.fp(mol)); hash.update('\n');
        });
      } else {
        const mol = prepared[i];
        hash.update(api.fp(mol)); hash.update('\n');
      }
    }
    const elapsed_ms = Number(process.hrtime.bigint() - start) / 1e6;
    return { elapsed_ms, parsed_rows: values.length,
      morgan_rows_sha256: lane === 'parse_only' ? null : hash.digest('hex') };
  } finally {
    for (const mol of prepared) api.free(mol);
  }
}
async function child() {
  const name = arg('--package');
  const lane = arg('--lane');
  check(lane === 'preflight' || LANES.includes(lane), 'unknown lane');
  const corpus = rows();
  const api = await load(name);
  if (lane === 'preflight') {
    process.stdout.write(JSON.stringify({ name, lane, ...preflight(api, corpus.values) }) + '\n');
    return;
  }
  warm(api, corpus.values, lane);
  const measure = timed(api, corpus.values, lane);
  process.stdout.write(JSON.stringify({
    name, lane, ...measure, init_ms: api.init_ms,
    peak_rss_bytes: process.resourceUsage().maxRSS * 1024,
    rss_after_bytes: process.memoryUsage().rss,
  }) + '\n');
}
function bootstrap(logRatios) {
  let state = 0x9e3779b9;
  const sample = [];
  for (let repeat = 0; repeat < 10000; repeat++) {
    let mean = 0;
    for (let i = 0; i < logRatios.length; i++) {
      state ^= state << 13; state ^= state >>> 17; state ^= state << 5;
      mean += logRatios[(state >>> 0) % logRatios.length];
    }
    sample.push(Math.exp(mean / logRatios.length));
  }
  return [q(sample, 0.025), q(sample, 0.975)];
}
function runChild(base, name, lane) {
  const result = spawnSync(process.execPath, [SELF, '--child', ...base, '--package', name, '--lane', lane], {
    encoding: 'utf8', maxBuffer: 8 * 1024 * 1024,
  });
  if (result.status !== 0) throw new Error(`${name}/${lane}: ${result.stderr || result.stdout}`);
  return JSON.parse(result.stdout.trim());
}
function comparison(base, left, right, preflights, lane) {
  const blocks = [];
  for (let block = 0; block < BLOCKS; block++) {
    const order = block % 2 === 0 ? [left, right, right, left] : [right, left, left, right];
    const samples = order.map(name => runChild(base, name, lane));
    const sum = name => samples.filter(x => x.name === name).reduce((a, x) => a + x.elapsed_ms, 0) / 2;
    const leftMs = sum(left), rightMs = sum(right);
    check(leftMs > 0 && rightMs > 0, 'nonpositive timing');
    for (const sample of samples) {
      check(sample.parsed_rows === preflights.compared_rows, 'timed row count mismatch');
      if (lane !== 'parse_only')
        check(sample.morgan_rows_sha256 === preflights.morgan_rows_sha256, 'timed Morgan output mismatch');
    }
    blocks.push({ index: block, order, samples, left_ms: leftMs, right_ms: rightMs,
      speed_ratio_left_over_right: rightMs / leftMs });
    if ((block + 1) % 5 === 0) process.stderr.write(`${left}/${right} ${lane}: ${block + 1}/20 blocks\n`);
  }
  const logs = blocks.map(x => Math.log(x.speed_ratio_left_over_right));
  const ci = bootstrap(logs);
  return { left, right, lane, blocks, summary: {
    median_speed_ratio_left_over_right: q(blocks.map(x => x.speed_ratio_left_over_right), 0.5),
    geometric_mean_ratio: Math.exp(logs.reduce((a, b) => a + b, 0) / logs.length),
    bootstrap_95_ci: ci,
    statistical_speed_signal: ci[0] > 1,
    practical_10pct_speed_win: ci[0] > 1.10,
    interval_crosses_parity: ci[0] <= 1 && ci[1] >= 1,
    left_init_median_ms: q(blocks.flatMap(x => x.samples.filter(s => s.name === left).map(s => s.init_ms)), 0.5),
    right_init_median_ms: q(blocks.flatMap(x => x.samples.filter(s => s.name === right).map(s => s.init_ms)), 0.5),
    left_peak_rss_median_bytes: q(blocks.flatMap(x => x.samples.filter(s => s.name === left).map(s => s.peak_rss_bytes)), 0.5),
    right_peak_rss_median_bytes: q(blocks.flatMap(x => x.samples.filter(s => s.name === right).map(s => s.peak_rss_bytes)), 0.5),
  } };
}
async function parent() {
  const output = resolve(arg('--output'));
  const corpus = rows();
  const names = ['v29', 'v30', 'rdkit'];
  const artifacts = {};
  for (const name of names) {
    const info = optionsFor(name);
    const paths = artifactPaths(name);
    const tarball = resolve(arg(`--${name}-tarball`));
    check(fileSha(tarball) === info.tarball, `${name} tarball hash mismatch`);
    artifacts[name] = { version: info.version, tarball_sha256: info.tarball,
      js_sha256: info.js, wasm_sha256: info.wasm };
  }
  const base = ['--corpus', resolve(arg('--corpus')), '--rows', String(corpus.count),
    ...names.flatMap(name => [`--${name}-dir`, resolve(arg(`--${name}-dir`))])];
  const preflightOutput = Object.fromEntries(names.map(name => [name, runChild(base, name, 'preflight')]));
  const reference = preflightOutput.v30.rows;
  for (const name of names) {
    check(preflightOutput[name].rows.length === corpus.count, `${name} preflight count mismatch`);
    check(preflightOutput[name].morgan_rows_sha256 === preflightOutput.v30.morgan_rows_sha256,
      `${name} preflight Morgan digest mismatch`);
    for (let i = 0; i < corpus.count; i++) {
      check(JSON.stringify(preflightOutput[name].rows[i]) === JSON.stringify(reference[i]),
        `${name} atom count/Morgan mismatch at row ${i}`);
    }
  }
  const atomCountSum = reference.reduce((a, row) => a + row[0], 0);
  const preflights = { compared_rows: corpus.count, exact_atom_and_morgan_rows: corpus.count,
    atom_count_sum: atomCountSum, morgan_rows_sha256: preflightOutput.v30.morgan_rows_sha256,
    by_artifact: Object.fromEntries(names.map(name => [name, preflightOutput[name].rows])) };
  const comparisons = [];
  for (const lane of LANES) {
    comparisons.push(comparison(base, 'v29', 'v30', preflights, lane));
    comparisons.push(comparison(base, 'v30', 'rdkit', preflights, lane));
  }
  const os = await import('node:os');
  const record = {
    schema: 'published-wasm-paired-speed/v1',
    runner_sha256: fileSha(SELF),
    environment: { node: process.version, platform: process.platform, arch: process.arch,
      os_release: os.release(), cpu_count: os.cpus().length,
      cpus: os.cpus().map(x => x.model).filter((x, i, a) => a.indexOf(x) === i) },
    corpus: { path: 'scripts/descriptor_census_corpus.smi', sha256: corpus.hash, rows: corpus.count },
    artifacts, preflight: preflights,
    protocol: { paired_blocks: BLOCKS, processes_per_block: 4,
      order: 'ABBA then BAAB alternating', warmup_rows: Math.min(20, corpus.count),
      confidence_interval: '10,000 deterministic bootstrap resamples of 20 paired log-ratios',
      practical_win_rule: 'lower 95% interval > 1.10; smaller effects remain diagnostic',
      lanes: LANES, boundary: 'one fresh Node process per artifact and observation',
      memory: 'process-wide peak RSS includes Node/WASM initialization and operations; not matched library allocations' },
    comparisons,
    interpretation: 'Node-only published-artifact measurement. No browser, cross-host, canonical-string, or library-allocation claim.',
  };
  writeFileSync(output, JSON.stringify(record, null, 2) + '\n');
  process.stdout.write(JSON.stringify({ output, summaries: comparisons.map(x => ({
    left: x.left, right: x.right, lane: x.lane, ...x.summary,
  })) }) + '\n');
}

if (process.argv.includes('--child')) await child();
else await parent();
