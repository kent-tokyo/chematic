#!/usr/bin/env python3
"""Run BioTransformer's public reaction rules through RDKit and chematic.

#734 asked for the roughly 1,200 BioTransformer metabolism rules to be run on
the same reactants in both engines. The public BioTransformer repository
(Wishartlab-openscience/Biotransformer) keeps its rule tables as JSON; this
script reads them from a local checkout or download (they are LGPL-3.0 with a
non-commercial clause and the environmental tables are CC BY-NC-SA, so they
are not vendored here; ``--rules`` files are pinned by SHA-256 in the output).

For each rule and reactant the script applies the SMIRKS with RDKit
``RunReactants`` and with ``chematic.run_smirks_checked(..., rdkit_compat=True)``
and compares the product sets after RDKit's own sanitize (an RDKit product
that fails ``SanitizeMol`` is not a product; this is the "after sanitize"
column of the #734 harness). Reactants are run as parsed and with explicit
hydrogens (``AddHs``), since most rules spell ``[H]`` atoms. Each reactant
is written as RDKit canonical SMILES and read again before either engine sees
it, so both start from the same molecule.

Outcome per (rule, reactant, hydrogen mode):

* ``both_none``: neither engine returns a sanitizable product;
* ``both_none_refused``: chematic refuses and RDKit has no product;
* ``exact``: the same set of product sets;
* ``chematic_refused``: typed refusal / unsupported / error from chematic
  where RDKit has products;
* ``rdkit_raw_unsanitizable``: RDKit returns raw products that all fail its
  own sanitize, chematic returns products;
* ``rdkit_truncated_subset``: RDKit stopped at ``--max-products`` and its
  set is a subset of chematic's;
* ``rdkit_resanitize_fails_same_products``: the sets differ only because
  RDKit's ``RemoveHs`` re-sanitizes a product that RDKit's own ``SanitizeMol``
  accepted and fails (an aromatic ring perceived across order-less ``~``
  bonds that its kekulizer cannot redo), or because RDKit writes the same
  hydrogen count with a different implicit-H flag: compared after a single
  sanitize with every atom's hydrogen count made explicit, the two sets agree;
* ``differ``: anything else, with both sets recorded.

Rules RDKit cannot parse, and rules with more than one reactant template,
are counted and skipped.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from collections import Counter
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")


def load_rules(paths: list[Path]) -> tuple[list[dict], list[dict]]:
    import json5

    rules, sources = [], []
    for path in paths:
        raw = path.read_bytes()
        sources.append({"path": path.name, "sha256": hashlib.sha256(raw).hexdigest()})
        text = raw.decode("utf-8")
        table = json5.loads(text[text.index("{"):], object_pairs_hook=lambda pairs: pairs)
        top = dict(table)
        entries = top.get("reactions", table)
        for name, body in entries:
            body = dict(body)
            if body.get("smirks"):
                rules.append({
                    "table": path.name,
                    "name": name,
                    "id": body.get("btmrID") or name,
                    "smirks": body["smirks"],
                })
    return rules, sources


def sanitized_key(mols) -> tuple[str, ...] | None:
    """Sorted canonical SMILES of the sanitized fragments, or None."""
    out = []
    for mol in mols:
        mol = Chem.Mol(mol)
        try:
            Chem.SanitizeMol(mol)
            mol = Chem.RemoveHs(mol)
        except Exception:
            return None
        smi = Chem.MolToSmiles(mol)
        out.extend(smi.split("."))
    return tuple(sorted(out))


def single_sanitize_key(mols) -> tuple[str, ...] | None:
    """Sorted canonical SMILES after one ``SanitizeMol`` (``RemoveHs`` without
    re-sanitizing), each atom's hydrogen count written explicitly."""
    out = []
    for mol in mols:
        mol = Chem.Mol(mol)
        try:
            Chem.SanitizeMol(mol)
            mol = Chem.RemoveHs(mol, sanitize=False)
        except Exception:
            return None
        for atom in mol.GetAtoms():
            h = atom.GetTotalNumHs()
            atom.SetNoImplicit(True)
            atom.SetNumExplicitHs(h)
        out.extend(Chem.MolToSmiles(mol).split("."))
    return tuple(sorted(out))


def rdkit_set(rxn, mol, max_products: int) -> tuple[set, int, list]:
    """Sanitized product sets, the raw product-set count and the raw sets."""
    values = set()
    raw = rxn.RunReactants((mol,), max_products)
    for product_set in raw:
        key = sanitized_key(product_set)
        if key is not None:
            values.add(key)
    return values, len(raw), raw


def to_rdkit(product):
    """An RDKit copy of a chematic product. RDKit's SMILES reader drops `~`
    on a ring closure, so order-less (zero-order) bonds are set from the MOL
    block's type-8 bonds, mapped through the SMILES atom order."""
    if "~" not in product.smiles:
        return Chem.MolFromSmiles(product.smiles, sanitize=False)
    smi, order = product.smiles_with_atom_order()
    mol = Chem.MolFromSmiles(smi, sanitize=False)
    if mol is None:
        return None
    rd_of = {chem: rd for rd, chem in enumerate(order)}
    lines = product.to_mol_block().splitlines()
    natoms, nbonds = int(lines[3][0:3]), int(lines[3][3:6])
    rw = Chem.RWMol(mol)
    for line in lines[4 + natoms: 4 + natoms + nbonds]:
        a, b, kind = int(line[0:3]) - 1, int(line[3:6]) - 1, int(line[6:9])
        if kind == 8:
            bond = rw.GetBondBetweenAtoms(rd_of[a], rd_of[b])
            if bond is not None:
                bond.SetBondType(Chem.BondType.UNSPECIFIED)
    return rw.GetMol()


def chematic_set(chematic, smirks: str, smiles: str) -> tuple[str, tuple | None, str | None]:
    try:
        reactant = chematic.from_smiles(smiles)
    except ValueError as exc:
        return "reactant_parse", None, str(exc)
    try:
        checked = chematic.run_smirks_checked(smirks, [reactant], rdkit_compat=True)
    except ValueError as exc:
        return "error", None, str(exc)
    status = checked["status"]
    if status in {"typed_refusal", "typed_unsupported"}:
        return status, None, checked.get("reason")
    values = set()
    single = set()
    for product_set in checked.get("products", []):
        mols = [to_rdkit(p) for p in product_set]
        if any(m is None for m in mols):
            continue
        key = sanitized_key(mols)
        if key is not None:
            values.add(key)
        key = single_sanitize_key(mols)
        if key is not None:
            single.add(key)
    return status, (values, single), None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rules", type=Path, nargs="+", required=True)
    ap.add_argument("--reactants", type=Path, required=True,
                    help="SMILES file (first column)")
    ap.add_argument("--every", type=int, default=1, help="take every n-th reactant row")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--max-products", type=int, default=1000)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--rows", type=Path, required=True, help="JSONL of non-exact rows")
    args = ap.parse_args()

    import chematic

    rules, sources = load_rules(args.rules)
    lines = [l.split()[0] for l in args.reactants.read_text().splitlines() if l.strip()]
    lines = lines[:: args.every]
    if args.limit:
        lines = lines[: args.limit]
    reactants = []
    for smi in lines:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        # Both engines read the same canonical SMILES. RDKit keeps the H
        # count of a bracket atom in the input (`[C]` written for a carbon
        # with four bonds) as fixed, which its canonical SMILES drops; a
        # template that removes a bond then gives a radical in RDKit only.
        smi = Chem.MolToSmiles(mol)
        mol = Chem.MolFromSmiles(smi)
        reactants.append((smi, mol, Chem.AddHs(mol)))

    counts = Counter()
    per_rule = {}
    started = time.time()
    with args.rows.open("w") as rows_out:
        for rule in rules:
            try:
                rxn = AllChem.ReactionFromSmarts(rule["smirks"])
            except ValueError:
                rxn = None
            if rxn is None:
                counts["rdkit_parse_error_rules"] += 1
                per_rule[rule["id"]] = {"rdkit_parse_error": True}
                continue
            if rxn.GetNumReactantTemplates() != 1:
                counts["multi_reactant_rules_skipped"] += 1
                per_rule[rule["id"]] = {"skipped": "reactant templates != 1"}
                continue
            rxn.Initialize()
            stats = Counter()
            for smi, mol, molh in reactants:
                for mode, rmol in (("implicit", mol), ("explicit_h", molh)):
                    want, raw_count, raw = rdkit_set(rxn, rmol, args.max_products)
                    status, got_pair, detail = chematic_set(chematic, rule["smirks"], Chem.MolToSmiles(rmol))
                    got, got_single = got_pair if got_pair is not None else (None, None)
                    if got is None:
                        outcome = "chematic_refused" if want else "both_none_refused"
                    elif not want and not got:
                        outcome = "both_none"
                    elif want == got:
                        outcome = "exact"
                    elif raw_count < args.max_products and got_single == (
                        {k for k in map(single_sanitize_key, raw) if k is not None}
                    ):
                        outcome = "rdkit_resanitize_fails_same_products"
                    elif raw_count >= args.max_products and want <= got:
                        # RDKit stopped at max_products (H-atom mapping
                        # permutations), so its set is a subset.
                        outcome = "rdkit_truncated_subset"
                    elif raw_count and not want:
                        # Every raw RDKit product fails RDKit's own sanitize.
                        outcome = "rdkit_raw_unsanitizable"
                    else:
                        outcome = "differ"
                    counts[f"{mode}:{outcome}"] += 1
                    stats[f"{mode}:{outcome}"] += 1
                    if outcome in {"differ", "chematic_refused", "rdkit_raw_unsanitizable",
                                   "rdkit_truncated_subset",
                                   "rdkit_resanitize_fails_same_products"}:
                        rows_out.write(json.dumps({
                            "rule": rule["id"], "name": rule["name"], "table": rule["table"],
                            "smirks": rule["smirks"], "reactant": smi, "mode": mode,
                            "outcome": outcome, "chematic_status": status, "detail": detail,
                            "rdkit": sorted(".".join(k) for k in want),
                            "chematic": sorted(".".join(k) for k in (got or [])),
                        }) + "\n")
            per_rule[rule["id"]] = dict(stats)
            counts["rules_run"] += 1
    summary = {
        "schema": "biotransformer-rule-corpus/v1",
        "rules_sources": sources,
        "rule_count": len(rules),
        "reactants": {"file": args.reactants.name,
                      "sha256": hashlib.sha256(args.reactants.read_bytes()).hexdigest(),
                      "every": args.every, "limit": args.limit, "count": len(reactants)},
        "rdkit": rdBase.rdkitVersion,
        "chematic": getattr(chematic, "__version__", None),
        "elapsed_seconds": round(time.time() - started, 1),
        "counts": dict(sorted(counts.items())),
        "per_rule": per_rule,
    }
    args.output.write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary["counts"], indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
