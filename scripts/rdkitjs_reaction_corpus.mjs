#!/usr/bin/env node
/**
 * Run SMIRKS rules over reactants with an RDKit.js (MinimalLib) build and
 * write each (rule, reactant, hydrogen mode)'s sanitized product sets.
 *
 * The RDKit Python 2026.09.1 wheel is not published, while RDKit.js
 * 2026.09.1 is: running this script once with RDKit.js 2026.03.6 and once
 * with 2026.09.1 on the same inputs shows which reaction outputs changed
 * between the two releases, the input to re-baselining chematic's RDKit
 * 2026.03.6 reaction profile.
 *
 * A raw product is sanitized by reading its MOL block back (RDKit.js exposes
 * no sanitize call), hydrogens removed, written as canonical SMILES; a product
 * set with a product RDKit cannot read back is dropped, like a product that
 * fails sanitize in the Python harness.
 *
 * Usage: node rdkitjs_reaction_corpus.mjs --rdkit DIST_DIR --rules RULES.json
 *          --reactants REACTANTS.json --output OUT.jsonl [--max-products 1000]
 *          [--rule-index I] [--start R] [--skip r:mode,...] [--progress FILE]
 *   RULES.json: [{"id", "smirks"}], REACTANTS.json: [smiles]
 *
 * RDKit.js can abort the WebAssembly instance on some inputs; the driver
 * (``rdkitjs_reaction_corpus_driver.py``) runs one rule per process
 * (``--rule-index``), reads ``--progress`` (the row being run) after an abort
 * and reruns the rule with that row in ``--skip``.
 */

import fs from 'node:fs';
import path from 'node:path';

function arg(name, fallback) {
  const i = process.argv.indexOf(name);
  if (i < 0) {
    if (fallback !== undefined) return fallback;
    throw new Error(`missing ${name}`);
  }
  return process.argv[i + 1];
}

const init = (await import(path.resolve(arg('--rdkit'), 'RDKit_minimal.js'))).default;
const RD = await init();
const allRules = JSON.parse(fs.readFileSync(arg('--rules')));
const ruleIndex = arg('--rule-index', null);
const rules = ruleIndex === null ? allRules : [allRules[Number(ruleIndex)]];
const skip = new Set((arg('--skip', '') || '').split(',').filter(Boolean));
const startRow = Number(arg('--start', '0'));
const progressPath = arg('--progress', null);
const progressFd = progressPath ? fs.openSync(progressPath, 'w') : null;
const reactants = JSON.parse(fs.readFileSync(arg('--reactants')));
const maxProducts = Number(arg('--max-products', '1000'));
// Synchronous writes: the run never yields to the event loop.
const fd = fs.openSync(arg('--output'), 'w');
const out = {
  write: line => fs.writeSync(fd, line),
  end: line => {
    fs.writeSync(fd, line);
    fs.closeSync(fd);
  },
};
out.write(JSON.stringify({ rdkit: RD.version(), rules: rules.length, reactants: reactants.length }) + '\n');

const molsFor = reactants.map(smi => {
  const plain = RD.get_mol(smi);
  const withH = RD.get_mol(smi);
  withH.add_hs_in_place();
  return { implicit: plain, explicit_h: withH };
});

function productSet(list) {
  const smiles = [];
  for (let i = 0; i < list.size(); i++) {
    const raw = list.at(i);
    let block = null;
    try {
      block = raw.get_molblock();
    } catch (e) {
      block = null;
    }
    raw.delete();
    if (!block) return null;
    const mol = RD.get_mol(block, JSON.stringify({ removeHs: true }));
    if (!mol || !mol.is_valid()) {
      if (mol) mol.delete();
      return null;
    }
    smiles.push(...mol.get_smiles().split('.'));
    mol.delete();
  }
  return smiles.sort().join('.');
}

const started = Date.now();
for (const rule of rules) {
  let rxn = null;
  try {
    rxn = RD.get_rxn(rule.smirks);
  } catch (e) {
    rxn = null;
  }
  if (!rxn) {
    out.write(JSON.stringify({ rule: rule.id, error: 'parse' }) + '\n');
    continue;
  }
  if (startRow === 0) out.write(JSON.stringify({ rule: rule.id, start: true }) + '\n');
  for (let r = startRow; r < reactants.length; r++) {
    for (const mode of ['implicit', 'explicit_h']) {
      if (skip.has(`${r}:${mode}`)) {
        out.write(JSON.stringify({ rule: rule.id, r, mode, error: 'abort' }) + '\n');
        continue;
      }
      if (progressFd !== null) {
        const mark = Buffer.from(`${r}:${mode}`.padEnd(32));
        fs.writeSync(progressFd, mark, 0, mark.length, 0);
      }
      const list = new RD.MolList();
      list.append(molsFor[r][mode]);
      let sets;
      let raw = 0;
      try {
        const products = rxn.run_reactants(list, maxProducts);
        raw = products.size();
        sets = new Set();
        for (let k = 0; k < raw; k++) {
          const pl = products.get(k);
          const key = productSet(pl);
          pl.delete();
          if (key !== null) sets.add(key);
        }
        products.delete();
      } catch (e) {
        sets = null;
      }
      list.delete();
      if (sets === null) {
        out.write(JSON.stringify({ rule: rule.id, r, mode, error: 'run' }) + '\n');
      } else if (raw > 0) {
        out.write(JSON.stringify({ rule: rule.id, r, mode, raw, sets: [...sets].sort() }) + '\n');
      }
    }
  }
  rxn.delete();
}
out.end(JSON.stringify({ seconds: (Date.now() - started) / 1000 }) + '\n');
