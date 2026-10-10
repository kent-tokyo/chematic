#!/usr/bin/env python3
"""Property-based SMIRKS fuzzing of chematic against RDKit (#734, #754).

The BioTransformer corpus runs about a thousand rules that people wrote. This
lane writes rules: from a reactant molecule it takes a connected set of one
to four atoms, writes a reactant template that matches them (element, H
count, degree, connectivity, charge, ring membership and aromaticity
primitives, joined with ``;``/``&``/``,``/``!``, and the bond primitives
``-``/``=``/``#``/``:``/``~``/``@``/``!@``), sometimes perturbs one primitive
so the template no longer matches, and applies one or two edits in the
product template:

* change a bond order, or break a bond;
* change a charge or an explicit H count;
* change an element;
* delete a mapped atom, or attach a new unmapped atom.

Each template is applied with RDKit ``RunReactants`` and with
``chematic.run_smirks_checked(..., rdkit_compat=True)`` to the molecule it
came from and to other reactants, implicit-H and with ``AddHs``. Outcomes
are classified exactly as in ``biotransformer_rule_corpus.py`` (products
compared after RDKit's own sanitize). Every ``differ`` and
``chematic_refused`` case is shrunk: primitives, edits and template atoms are
dropped, and the smallest reactant from the pool is kept, while the outcome
stays the same, so a reported probe is minimal.

The generator is seeded; the same ``--seed`` and inputs give the same
templates, so a run is reproducible and its probes can be replayed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import random
import sys
import time
from collections import Counter
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

sys.path.insert(0, str(Path(__file__).resolve().parent))
from biotransformer_rule_corpus import InProcess, Worker, chematic_set, classify, rdkit_set, sanitized_key  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REACTANTS = ROOT / "validation/biotransformer-requests-400.tsv"
REPORTED = {"differ", "chematic_refused"}


# ---------------------------------------------------------------------------
# Template generation
# ---------------------------------------------------------------------------

def pick_atoms(mol: Chem.Mol, rng: random.Random, size: int) -> list[int]:
    start = rng.randrange(mol.GetNumAtoms())
    chosen = [start]
    while len(chosen) < size:
        frontier = sorted({n.GetIdx() for i in chosen for n in mol.GetAtomWithIdx(i).GetNeighbors()} - set(chosen))
        if not frontier:
            break
        chosen.append(rng.choice(frontier))
    return chosen


def atom_primitives(atom: Chem.Atom, rng: random.Random) -> list[str]:
    """True primitives for this atom, the element first."""
    z, sym = atom.GetAtomicNum(), atom.GetSymbol()
    element = rng.choice([f"#{z}", sym.lower() if atom.GetIsAromatic() and sym in "CNOPS" else sym])
    if element in ("C", "N", "O", "S", "P") and atom.GetIsAromatic():
        element = f"#{z}"
    prims = [element]
    pool = [
        f"H{atom.GetTotalNumHs()}",
        f"D{atom.GetDegree()}",
        f"X{atom.GetTotalDegree()}",
        "R" if atom.IsInRing() else "!R",
        "a" if atom.GetIsAromatic() else "A",
        f"+{atom.GetFormalCharge()}" if atom.GetFormalCharge() >= 0 else f"-{-atom.GetFormalCharge()}",
    ]
    if atom.IsInRing():
        pool.append(f"r{min(s for s in range(3, 9) if atom.IsInRingSize(s))}" if any(atom.IsInRingSize(s) for s in range(3, 9)) else "R")
    prims += rng.sample(pool, k=rng.randint(0, 3))
    return prims


def false_primitive(atom: Chem.Atom, rng: random.Random) -> str:
    """A primitive this atom does not satisfy."""
    options = [
        f"H{atom.GetTotalNumHs() + 1}",
        f"D{atom.GetDegree() + 1}",
        "!R" if atom.IsInRing() else "R",
        "A" if atom.GetIsAromatic() else "a",
        f"#{atom.GetAtomicNum() + 1}",
    ]
    return rng.choice(options)


def join_primitives(prims: list[str], rng: random.Random) -> str:
    if len(prims) == 1:
        return prims[0]
    out = prims[0]
    for p in prims[1:]:
        out += rng.choice([";", "&", ""]) + p
    if rng.random() < 0.15:
        # An alternative that the atom also satisfies (an OR with itself).
        out = f"{out},{prims[0]}" if ";" not in out and "&" not in out else out
    return out


def bond_primitive(bond: Chem.Bond, rng: random.Random) -> str:
    order = {Chem.BondType.SINGLE: "-", Chem.BondType.DOUBLE: "=", Chem.BondType.TRIPLE: "#",
             Chem.BondType.AROMATIC: ":"}.get(bond.GetBondType(), "~")
    ring = "@" if bond.IsInRing() else "!@"
    return rng.choice([order, order, "~", f"{order};{ring}", ring])


def write_graph(atoms: list[int], labels: dict[int, str], bonds: dict[tuple[int, int], str]) -> str:
    """A SMARTS/SMILES string for the given atoms (connected), ring closures
    for bonds that close a cycle."""
    adj: dict[int, list[int]] = {a: [] for a in atoms}
    for (a, b) in bonds:
        adj[a].append(b)
        adj[b].append(a)
    seen: set[int] = set()
    closures: dict[tuple[int, int], int] = {}
    tree: set[tuple[int, int]] = set()

    def visit(a: int, parent: int | None) -> None:
        seen.add(a)
        for b in adj[a]:
            if b == parent:
                continue
            key = (min(a, b), max(a, b))
            if b in seen:
                if key not in tree and key not in closures:
                    closures[key] = len(closures) + 1
            else:
                tree.add(key)
                visit(b, a)

    order: list[int] = []
    components = []
    for a in atoms:
        if a not in seen:
            visit(a, None)
            components.append(a)
    seen.clear()

    def emit(a: int, parent: int | None) -> str:
        seen.add(a)
        out = labels[a]
        for key, digit in closures.items():
            if a in key:
                other = key[0] if key[1] == a else key[1]
                bond = bonds[key] if other not in seen else ""
                out += f"{bond}%{digit + 9}" if digit + 9 > 9 else f"{bond}{digit}"
        children = [b for b in adj[a] if b != parent and (min(a, b), max(a, b)) in tree and b not in seen]
        for i, b in enumerate(children):
            branch = bonds[(min(a, b), max(a, b))] + emit(b, a)
            out += f"({branch})" if i < len(children) - 1 else branch
        return out

    return ".".join(emit(c, None) for c in components)


def generate(mol: Chem.Mol, rng: random.Random) -> dict | None:
    size = rng.randint(1, 4)
    atoms = pick_atoms(mol, rng, size)
    maps = {a: i + 1 for i, a in enumerate(atoms)}
    prims = {a: atom_primitives(mol.GetAtomWithIdx(a), rng) for a in atoms}
    if rng.random() < 0.15:
        a = rng.choice(atoms)
        prims[a].append(false_primitive(mol.GetAtomWithIdx(a), rng))
    bonds = {}
    for i, a in enumerate(atoms):
        for b in atoms[i + 1:]:
            bond = mol.GetBondBetweenAtoms(a, b)
            if bond is not None:
                bonds[(min(a, b), max(a, b))] = bond_primitive(bond, rng)

    product_labels = {a: f"[*:{maps[a]}]" for a in atoms}
    product_bonds = {}
    for key in bonds:
        bond = mol.GetBondBetweenAtoms(*key)
        product_bonds[key] = {Chem.BondType.SINGLE: "-", Chem.BondType.DOUBLE: "=",
                              Chem.BondType.TRIPLE: "#", Chem.BondType.AROMATIC: ":"}.get(bond.GetBondType(), "~")
    edits = []
    extra = []
    for _ in range(rng.randint(1, 2)):
        kind = rng.choice(["bond_order", "break", "charge", "hcount", "element", "delete", "add"])
        if kind in ("bond_order", "break") and product_bonds:
            key = rng.choice(sorted(product_bonds))
            if kind == "break":
                del product_bonds[key]
            else:
                product_bonds[key] = {"-": "=", "=": "-", "#": "=", ":": "-", "~": "-"}[product_bonds[key]]
        elif kind == "charge":
            a = rng.choice([a for a in atoms if a in product_labels])
            sym = mol.GetAtomWithIdx(a).GetSymbol()
            charge = rng.choice(["+", "-", "+0"])
            product_labels[a] = f"[{sym}{charge}:{maps[a]}]"
        elif kind == "hcount":
            a = rng.choice([a for a in atoms if a in product_labels])
            sym = mol.GetAtomWithIdx(a).GetSymbol()
            product_labels[a] = f"[{sym}H{rng.randint(0, 3)}:{maps[a]}]"
        elif kind == "element":
            a = rng.choice([a for a in atoms if a in product_labels])
            product_labels[a] = f"[{rng.choice(['C', 'N', 'O', 'S'])}:{maps[a]}]"
        elif kind == "delete" and len([a for a in atoms if a in product_labels]) > 1:
            a = rng.choice([a for a in atoms[1:] if a in product_labels] or [None])
            if a is None:
                continue
            product_labels.pop(a)
            product_bonds = {k: v for k, v in product_bonds.items() if a not in k}
        elif kind == "add":
            extra.append((rng.choice(atoms), rng.choice(["O", "C", "N", "Cl"])))
        else:
            continue
        edits.append(kind)
    kept = [a for a in atoms if a in product_labels]
    if not kept:
        return None
    # New atoms take indices past the molecule's.
    next_index = mol.GetNumAtoms()
    for anchor, element in extra:
        if anchor not in product_labels:
            continue
        product_labels[next_index] = element
        product_bonds[(anchor, next_index)] = "-"
        kept.append(next_index)
        next_index += 1
    reactant_labels = {a: f"[{join_primitives(prims[a], rng)}:{maps[a]}]" for a in atoms}
    try:
        reactant = write_graph(atoms, reactant_labels, bonds)
        product = write_graph(kept, product_labels, product_bonds)
    except KeyError:
        return None
    smirks = f"{reactant}>>{product}"
    if "." in reactant:
        return None
    try:
        rxn = AllChem.ReactionFromSmarts(smirks)
    except ValueError:
        return None
    if rxn is None or rxn.GetNumReactantTemplates() != 1:
        return None
    return {"smirks": smirks, "edits": edits}


# ---------------------------------------------------------------------------
# Running and shrinking
# ---------------------------------------------------------------------------

def raw_key(mols) -> tuple[str, ...]:
    """Sorted canonical SMILES of the unsanitized fragments (aromatic flags
    and hydrogens as the engine left them, explicit H atoms removed)."""
    out = []
    for mol in mols:
        mol = Chem.Mol(mol)
        mol.UpdatePropertyCache(strict=False)
        try:
            mol = Chem.RemoveHs(mol, sanitize=False)
        except Exception:
            pass
        out.extend(Chem.MolToSmiles(mol).split("."))
    return tuple(sorted(out))


def outcome(backend, smirks: str, smiles: str, explicit_h: bool, max_products: int) -> tuple[str, dict]:
    """The ``classify`` outcome, refined for two decided cases:

    * ``implicit_h_equivalent``: an explicit-H reactant where chematic gives
      exactly RDKit's products for the implicit-H form (chematic's
      documented policy; RDKit keeps the H atoms as substituents);
    * ``same_graph_kekule_order``: the raw products are the same graphs and
      differ only after RDKit's kekulization, which depends on atom order
      (``[c+]1ccccc1C`` and ``Cc1[c+]cccc1`` kekulize to different
      structures in RDKit itself).
    """
    rxn = AllChem.ReactionFromSmarts(smirks)
    rxn.Initialize()
    implicit = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(implicit) if explicit_h else implicit
    want, raw_count, raw = rdkit_set(rxn, mol, max_products)
    status, got_pair, detail = chematic_set(backend, smirks, Chem.MolToSmiles(mol))
    result = classify(want, raw_count, raw, got_pair, max_products)
    got = got_pair[0] if got_pair is not None else None
    if result == "differ" and explicit_h and got is not None:
        want_implicit, _, _ = rdkit_set(rxn, implicit, max_products)
        if got == want_implicit:
            result = "implicit_h_equivalent"
    if result == "differ" and got is not None:
        response = backend.run({"smirks": smirks, "smiles": Chem.MolToSmiles(mol)})
        chematic_raw = set()
        for product_set in response["products"] or []:
            mols = [Chem.MolFromSmiles(p["smiles"], sanitize=False) for p in product_set]
            if all(m is not None for m in mols):
                chematic_raw.add(raw_key(mols))
        rdkit_raw = {raw_key(ps) for ps in raw if sanitized_key(ps) is not None}
        if chematic_raw and chematic_raw == rdkit_raw:
            result = "same_graph_kekule_order"
    return result, {
        "chematic_status": status, "detail": detail,
        "rdkit": sorted(".".join(k) for k in want),
        "chematic": sorted(".".join(k) for k in (got or [])),
    }


def explicit_h_residual_facets(result: str, edits: list[str], detail: dict) -> set[str]:
    """Describe, without reclassifying, why an explicit-H residual needs review.

    Facets intentionally overlap. A deletion can both leave a radical and fail
    the final valence/sanitize check; retaining both facts is more useful than
    forcing the case into one optimistic root cause.
    """
    if result not in REPORTED:
        return set()
    facets = {f"edit:{edit}" for edit in set(edits)}
    if "delete" in edits:
        facets.add("atom_deletion_after_explicit_h")
    if "break" in edits:
        facets.add("bond_break_after_explicit_h")
    if result == "chematic_refused":
        facets.add("resanitize_or_valence_refusal")

    for backend in ("rdkit", "chematic"):
        for product_set in detail.get(backend, []):
            for smiles in product_set.split("."):
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    mol = Chem.MolFromSmiles(smiles, sanitize=False)
                if mol is not None and any(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms()):
                    facets.add("radical_product")
                    return facets
    return facets


def simplifications(smirks: str) -> list[str]:
    """Templates with one primitive, bond primitive or bracket term fewer."""
    out = []
    for sep in (";", "&", ","):
        start = 0
        while (i := smirks.find(sep, start)) != -1:
            end = i + 1
            while end < len(smirks) and smirks[end] not in ";&,:]":
                end += 1
            if smirks[end - 1:end] != "" and smirks[end:end + 1] in (";", "&", ",", ":", "]"):
                out.append(smirks[:i] + smirks[end:])
            start = i + 1
    for a, b in (("-;!@", "-"), ("-;@", "-"), ("=;!@", "="), ("=;@", "="), (":;@", ":"), ("~", "")):
        if a in smirks:
            out.append(smirks.replace(a, b, 1))
    return [s for s in dict.fromkeys(out) if s != smirks]


def shrink(backend, case: dict, pool: list[str], max_products: int) -> dict:
    target = case["outcome"]
    smirks, smiles, explicit_h = case["smirks"], case["reactant"], case["mode"] == "explicit_h"

    def same(s: str, smi: str) -> bool:
        try:
            rxn = AllChem.ReactionFromSmarts(s)
        except ValueError:
            return False
        if rxn is None or rxn.GetNumReactantTemplates() != 1:
            return False
        try:
            return outcome(backend, s, smi, explicit_h, max_products)[0] == target
        except Exception:
            return False

    changed = True
    while changed:
        changed = False
        for candidate in simplifications(smirks):
            if same(candidate, smiles):
                smirks, changed = candidate, True
                break
    for smi in sorted(pool, key=lambda s: (Chem.MolFromSmiles(s).GetNumAtoms(), s))[:60]:
        if Chem.MolFromSmiles(smi).GetNumAtoms() >= Chem.MolFromSmiles(smiles).GetNumAtoms():
            break
        if same(smirks, smi):
            smiles = smi
            break
    _, detail = outcome(backend, smirks, smiles, explicit_h, max_products)
    return {"smirks": smirks, "reactant": smiles, "mode": case["mode"], "outcome": target, **detail}


def load_pool(path: Path) -> list[str]:
    pool = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        mol = Chem.MolFromSmiles(line.split("\t")[0].split()[0])
        if mol is not None:
            pool.append(Chem.MolToSmiles(mol))
    return list(dict.fromkeys(pool))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--reactants", type=Path, default=DEFAULT_REACTANTS,
                    help="SMILES (first column); default the 400 BioTransformer reactants")
    ap.add_argument("--templates", type=int, default=2000)
    ap.add_argument("--others", type=int, default=3, help="other reactants each template is applied to")
    ap.add_argument("--seed", type=int, default=754)
    ap.add_argument("--max-products", type=int, default=1000)
    ap.add_argument("--no-shrink", action="store_true")
    ap.add_argument("--chematic-python", default=None,
                    help="run chematic in this interpreter (a published wheel for a Python RDKit has no wheel for)")
    ap.add_argument("--expected", type=Path, default=None,
                    help="validation/smirks_property_fuzz_expected.json: fail when a count differs")
    ap.add_argument("--output", type=Path, required=True, help="summary JSON")
    ap.add_argument("--rows", type=Path, required=True, help="JSONL of reported (shrunk) cases")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    pool = load_pool(args.reactants)
    backend = Worker(args.chematic_python) if args.chematic_python else InProcess()
    counts: Counter = Counter()
    edits: Counter = Counter()
    explicit_h_residuals: Counter = Counter()
    reported: list[dict] = []
    probes: dict[tuple, dict] = {}
    started = time.time()
    generated = 0
    attempts = 0
    while generated < args.templates and attempts < args.templates * 20:
        attempts += 1
        source = rng.choice(pool)
        template = generate(Chem.MolFromSmiles(source), rng)
        if template is None:
            continue
        generated += 1
        edits.update(template["edits"])
        targets = [source] + rng.sample(pool, k=min(args.others, len(pool)))
        for smi in targets:
            for explicit_h in (False, True):
                result, detail = outcome(backend, template["smirks"], smi, explicit_h, args.max_products)
                mode = "explicit_h" if explicit_h else "implicit"
                counts[f"{mode}:{result}"] += 1
                if explicit_h:
                    explicit_h_residuals.update(
                        explicit_h_residual_facets(result, template["edits"], detail)
                    )
                if result in REPORTED:
                    case = {"smirks": template["smirks"], "reactant": smi, "mode": mode,
                            "outcome": result, "edits": template["edits"], **detail}
                    reported.append(case)
    rows = []
    for case in reported:
        small = case if args.no_shrink else shrink(backend, case, pool, args.max_products)
        key = (small["smirks"], small["reactant"], small["mode"])
        if key not in probes:
            probes[key] = small
            rows.append({"original": {k: case[k] for k in ("smirks", "reactant", "mode", "edits")}, **small})
    with args.rows.open("w", encoding="utf-8") as stream:
        for row in rows:
            stream.write(json.dumps(row) + "\n")
    summary = {
        "schema": "smirks-property-fuzz/v1",
        "seed": args.seed,
        "templates": generated,
        "attempts": attempts,
        "applications": sum(counts.values()),
        "reactants": {"file": args.reactants.name, "sha256": hashlib.sha256(args.reactants.read_bytes()).hexdigest(),
                      "count": len(pool)},
        "rdkit": rdBase.rdkitVersion,
        "chematic": backend.version,
        "elapsed_seconds": round(time.time() - started, 1),
        "edits": dict(sorted(edits.items())),
        "counts": dict(sorted(counts.items())),
        "explicit_h_residual_facets": dict(sorted(explicit_h_residuals.items())),
        "reported_cases": len(reported),
        "distinct_minimized_probes": len(rows),
    }
    if args.expected is not None:
        expected = json.loads(args.expected.read_text(encoding="utf-8"))
        same_run = (expected["seed"], expected["templates"], expected["rdkit"]) == (
            args.seed, generated, summary["rdkit"])
        if not same_run:
            raise SystemExit(f"expected counts are for seed {expected['seed']}, "
                             f"{expected['templates']} templates and RDKit {expected['rdkit']}")
        problems = [f"{k}: expected {expected['counts'].get(k, 0)}, got {counts.get(k, 0)}"
                    for k in sorted(set(expected["counts"]) | set(counts))
                    if expected["counts"].get(k, 0) != counts.get(k, 0)]
        expected_facets = expected.get("explicit_h_residual_facets")
        if expected_facets is not None:
            actual_facets = summary["explicit_h_residual_facets"]
            problems.extend(
                f"explicit_h_residual_facets.{key}: expected "
                f"{expected_facets.get(key, 0)}, got {actual_facets.get(key, 0)}"
                for key in sorted(set(expected_facets) | set(actual_facets))
                if expected_facets.get(key, 0) != actual_facets.get(key, 0)
            )
        summary["expected"] = str(args.expected)
        summary["differences"] = problems
    args.output.write_text(json.dumps(summary, indent=1) + "\n", encoding="utf-8")
    print(json.dumps({k: summary[k] for k in ("templates", "applications", "counts", "distinct_minimized_probes")}, indent=1))
    if summary.get("differences"):
        print("\n".join(summary["differences"]), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
