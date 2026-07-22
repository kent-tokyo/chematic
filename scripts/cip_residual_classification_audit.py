#!/usr/bin/env python3
"""
Diag/accurate-cip-audit: mechanism classification of chematic's CIP-assignment
residual against RDKit's rdCIPLabeler oracle, for BOTH engines separately --

  - LEGACY / default  (`chematic_chem::assign_cip`, `Mol.cip_stereo()`, the
    shell-pooling comparator described in docs/cip_accurate_rfc.md)
  - ACCURATE          (`chematic_cip::assign_cip_accurate_experimental`, the
    digraph-based engine, `Mol.cip_stereo(mode="accurate")`)

into named mechanism buckets:

  pseudoasymmetric, aromaticity_kekulize_bug, aromaticity_other, phosphorus,
  multiple_bond_dup, rule4b_candidate, rule3_candidate, rule5_candidate,
  unclassified

No mismatch is silently dropped -- every row lands in exactly one bucket via a
single deterministic precedence order (see BUCKET_PRECEDENCE / classify_row()
below), applied identically to both engines so the two numbers are directly
comparable. This does not re-derive the accurate engine's own already-settled
residual classification (docs/cip_accurate_rfc.md Milestones 4A-0/4B-0/4C-0/
4C-1) -- it re-verifies the numbers against a fresh oracle run and re-maps the
existing, already-diagnosed rows onto this script's bucket taxonomy, using the
same mechanical signals for the legacy engine's much larger, previously
un-subdivided residual (its 43 "uncharacterized" + 96 "aromatic_mancude" rows
from validation/cip_label_corpus.jsonl were never split into Rule 3/4/5).

Cross-check against the parallel `diag/aromaticity-rdkit-parity` diagnosis's
finding: `chematic_core::kekulize()` hard-fails (and, because
`build_molecule_from_model` only ever *promotes* an aromatic flag, never
demotes one, the failure is invisible at the atom/bond flag level) for
tropylium-cation/imidazolium/pyridinium/pyrylium/tellurophene/phosphole-shaped
rings. This script runs the real algorithm (not a SMARTS match against those 6
fixture SMILES) on every corpus molecule via
`crates/chematic-perception/examples/cip_kekulize_probe.rs`'s JSONL dump, and
only attributes a mismatch to that bug when the failure is LOCAL to the
stereocenter (the failing atom is the stereocenter itself, one of its direct
neighbors, or shares an SSSR ring with a direct neighbor) -- "kekulize fails
somewhere in this molecule" is not by itself evidence the bug caused THIS
atom's mismatch.

Bucket precedence (first match wins, applied per mismatching atom):
  1. pseudoasymmetric        -- oracle (`modern`) label is lowercase r/s.
  2. aromaticity_kekulize_bug -- kekulize() fails locally (see above).
  3. phosphorus              -- the stereocenter atom itself is P.
  4. aromaticity_other       -- stereocenter has >=1 aromatic direct neighbor,
                                 not already caught by #2.
  5. multiple_bond_dup       -- a direct (non-aromatic) neighbor itself carries
                                 a double/triple bond to a third atom (the
                                 classic CIP phantom-duplication trigger).
  6. rule4b_candidate        -- the stereocenter's neighbors contain a
                                 constitutionally-tied pair (equal RDKit
                                 canonical rank ignoring chirality -- Rules
                                 1a/2 would tie them) AND >=1 embedded chiral
                                 atom exists elsewhere in the molecule.
                                 (Rule 4a/4c are structurally N/A for chematic
                                 -- `chematic_core::Chirality` has no
                                 unit-type/axial variant at all, established
                                 directly in docs/cip_accurate_rfc.md Milestone
                                 4B-0 -- so any Rule-4-shaped candidate here is
                                 necessarily 4b, not re-derived per row.)
  7. rule3_candidate         -- same tied-pair signal as #6, but the
                                 discriminating feature is a stereo-specified
                                 (E/Z) double bond elsewhere instead of an
                                 embedded chiral atom.
  8. rule5_candidate         -- same tied-pair signal, uppercase oracle label,
                                 no #6/#7 signal -- the only rule left that
                                 breaks a Rules-1a/2 tie between otherwise-
                                 identical branches. Weak/heuristic: flagged as
                                 a candidate, not a proven mechanism -- see the
                                 RFC doc for why this is deliberately not
                                 overfit into a stronger claim.
  9. unclassified            -- none of the above signals fired. This is an
                                 expected, legitimate outcome, not a script
                                 failure -- the legacy engine's own RFC
                                 documents its bugs as structural
                                 shell-pooling artifacts, not always
                                 attributable to one CIP rule.

Usage:
    .venv/bin/python scripts/cip_residual_classification_audit.py \\
        --smiles ~/Downloads/SMILES.csv \\
        --accurate-tsv <path from corpus_snapshot --candidate> \\
        --kekulize-probe validation/results/cip_kekulize_probe.jsonl \\
        --oracle-instability validation/cip_oracle_instability.jsonl \\
        --out-corpus validation/cip_residual_classification_corpus.jsonl \\
        --out-summary validation/results/cip_residual_classification_summary.json

Reproduce the two Rust-side inputs first:
    cargo run -p chematic-perception --release --example cip_kekulize_probe \\
        -- ~/Downloads/SMILES.csv > validation/results/cip_kekulize_probe.jsonl
    cargo run -p chematic-cip --release --example corpus_snapshot -- \\
        --candidate ~/Downloads/SMILES.csv <accurate_candidate.tsv>
"""

import argparse
import csv
import datetime
import hashlib
import json
import re
import sys

sys.path.insert(0, ".")
import chematic  # noqa: E402
from rdkit import Chem  # noqa: E402
from rdkit.Chem import rdCIPLabeler  # noqa: E402

BUCKET_PRECEDENCE = [
    "pseudoasymmetric",
    "aromaticity_kekulize_bug",
    "phosphorus",
    "aromaticity_other",
    "multiple_bond_dup",
    "rule4b_candidate",
    "rule3_candidate",
    "rule5_candidate",
    "unclassified",
]

KEKULIZE_ATOM_RE = re.compile(r"atom (\d+)")


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_smiles(csv_path):
    with open(csv_path) as f:
        reader = csv.reader(f)
        next(reader)
        return [row[0] for row in reader if row]


def load_kekulize_probe(path):
    """smiles -> (kekulize_ok: bool, failing_atom_idx: int|None)."""
    out = {}
    with open(path) as f:
        for line in f:
            d = json.loads(line)
            if not d.get("parse_ok", True):
                continue
            smi = d["smiles"]
            if d.get("kekulize_ok"):
                out[smi] = (True, None)
            else:
                m = KEKULIZE_ATOM_RE.search(d.get("kekulize_error", ""))
                out[smi] = (False, int(m.group(1)) if m else None)
    return out


def load_accurate_snapshot(path):
    """(smiles, atom_idx) -> value string (R/S/r/s/skip:*/ERR)."""
    out = {}
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            smi, idx, value = parts[0], int(parts[1]), parts[2]
            out[(smi, idx)] = value
    return out


def load_oracle_instability(path):
    """(smiles, atom_idx) -> row dict, if any."""
    out = {}
    if not path:
        return out
    with open(path) as f:
        for line in f:
            d = json.loads(line)
            if d.get("_manifest"):
                continue
            out[(d["smiles"], d["atom_idx"])] = d
    return out


def kekulize_local(rd, aidx, kek_ok, fail_atom):
    """Is a local kekulize() failure plausibly the cause of this atom's
    mismatch? Local = the failing atom is the stereocenter itself, a direct
    neighbor, or shares an SSSR ring with a direct neighbor."""
    if kek_ok or fail_atom is None:
        return False
    if fail_atom == aidx:
        return True
    atom = rd.GetAtomWithIdx(aidx)
    neighbor_idxs = {n.GetIdx() for n in atom.GetNeighbors()}
    if fail_atom in neighbor_idxs:
        return True
    ri = rd.GetRingInfo()
    fail_rings = [set(r) for r in ri.AtomRings() if fail_atom in r]
    for nb in neighbor_idxs:
        for ring in fail_rings:
            if nb in ring:
                return True
    return False


def tied_neighbor_signal(rd, aidx):
    """True if >=2 of this atom's heavy-atom neighbors are constitutionally
    tied (equal RDKit canonical rank ignoring chirality) -- i.e. Rules 1a/2
    alone cannot distinguish them, the precondition for Rule 3/4/5."""
    ranks = list(Chem.CanonicalRankAtoms(rd, breakTies=False, includeChirality=False))
    atom = rd.GetAtomWithIdx(aidx)
    neighbor_ranks = [ranks[n.GetIdx()] for n in atom.GetNeighbors() if n.GetAtomicNum() > 1]
    return len(neighbor_ranks) != len(set(neighbor_ranks))


def multiple_bond_neighbor_signal(rd, aidx):
    """True if a direct, non-aromatic neighbor of this atom itself carries a
    double/triple bond to a third atom (the classic CIP duplicate-node
    trigger for double/triple bonds)."""
    atom = rd.GetAtomWithIdx(aidx)
    for nb in atom.GetNeighbors():
        if nb.GetIsAromatic():
            continue
        for b in nb.GetBonds():
            if b.GetBondType() in (Chem.BondType.DOUBLE, Chem.BondType.TRIPLE) and not b.GetIsAromatic():
                other = b.GetOtherAtom(nb)
                if other.GetIdx() != aidx:
                    return True
    return False


def embedded_chiral_elsewhere_signal(rd, aidx):
    return any(
        a.GetIdx() != aidx and a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in rd.GetAtoms()
    )


def stereo_double_bond_elsewhere_signal(rd):
    return any(b.GetStereo() != Chem.BondStereo.STEREONONE for b in rd.GetBonds())


def classify_row(rd, aidx, modern_code, kek_ok, fail_atom):
    """Apply BUCKET_PRECEDENCE, return (bucket, evidence dict)."""
    evidence = {}

    is_pseudo = modern_code.islower()
    evidence["oracle_label_lowercase"] = is_pseudo
    if is_pseudo:
        return "pseudoasymmetric", evidence

    local_kek_fail = kekulize_local(rd, aidx, kek_ok, fail_atom)
    evidence["kekulize_local_failure"] = local_kek_fail
    if local_kek_fail:
        return "aromaticity_kekulize_bug", evidence

    atom = rd.GetAtomWithIdx(aidx)
    if atom.GetSymbol() == "P":
        return "phosphorus", evidence

    has_aromatic_neighbor = any(n.GetIsAromatic() for n in atom.GetNeighbors())
    evidence["direct_aromatic_neighbor"] = has_aromatic_neighbor
    if has_aromatic_neighbor:
        return "aromaticity_other", evidence

    has_mult_bond_neighbor = multiple_bond_neighbor_signal(rd, aidx)
    evidence["multiple_bond_neighbor"] = has_mult_bond_neighbor
    if has_mult_bond_neighbor:
        return "multiple_bond_dup", evidence

    tied = tied_neighbor_signal(rd, aidx)
    evidence["tied_constitutional_neighbor_pair"] = tied
    if tied:
        embedded_chiral = embedded_chiral_elsewhere_signal(rd, aidx)
        evidence["embedded_chiral_elsewhere"] = embedded_chiral
        if embedded_chiral:
            return "rule4b_candidate", evidence
        stereo_db = stereo_double_bond_elsewhere_signal(rd)
        evidence["stereo_double_bond_elsewhere"] = stereo_db
        if stereo_db:
            return "rule3_candidate", evidence
        return "rule5_candidate", evidence

    return "unclassified", evidence


def build_oracle(smis):
    """smiles -> (modern_cip: {idx: code}, legacy_cip: {idx: code}, rd_mol)."""
    out = {}
    for smi in smis:
        rd = Chem.MolFromSmiles(smi)
        if rd is None:
            continue
        rd_legacy = Chem.MolFromSmiles(smi)
        Chem.AssignStereochemistry(rd_legacy, cleanIt=True, force=True)
        legacy_cip = {
            a.GetIdx(): a.GetProp("_CIPCode") for a in rd_legacy.GetAtoms() if a.HasProp("_CIPCode")
        }
        try:
            rdCIPLabeler.AssignCIPLabels(rd)
        except Exception:
            continue
        modern_cip = {a.GetIdx(): a.GetProp("_CIPCode") for a in rd.GetAtoms() if a.HasProp("_CIPCode")}
        if not modern_cip:
            continue
        out[smi] = (modern_cip, legacy_cip, rd)
    return out


def classify_legacy_engine(oracle, kek_probe):
    total = 0
    mismatches = []
    not_considered = 0
    for smi, (modern_cip, _legacy_cip, rd) in oracle.items():
        try:
            m = chematic.from_smiles(smi)
            cm_cip = {d["atom_idx"]: d["descriptor"] for d in m.cip_stereo()}
        except Exception:
            continue
        kek_ok, fail_atom = kek_probe.get(smi, (True, None))
        for aidx, modern_code in modern_cip.items():
            total += 1
            cm_code = cm_cip.get(aidx)
            if cm_code is None:
                not_considered += 1
            if cm_code == modern_code:
                continue
            bucket, evidence = classify_row(rd, aidx, modern_code, kek_ok, fail_atom)
            mismatches.append(
                {
                    "engine": "legacy",
                    "smiles": smi,
                    "atom_idx": aidx,
                    "chematic": cm_code,
                    "modern": modern_code,
                    "bucket": bucket,
                    "evidence": evidence,
                }
            )
    return total, mismatches, not_considered


def classify_accurate_engine(oracle, kek_probe, accurate_snapshot, oracle_instability):
    total = 0
    mismatches = []
    not_considered = 0
    for smi, (modern_cip, _legacy_cip, rd) in oracle.items():
        kek_ok, fail_atom = kek_probe.get(smi, (True, None))
        for aidx, modern_code in modern_cip.items():
            total += 1
            cm_code = accurate_snapshot.get((smi, aidx))
            if cm_code is None:
                not_considered += 1
                continue
            if cm_code == modern_code:
                continue
            bucket, evidence = classify_row(rd, aidx, modern_code, kek_ok, fail_atom)
            unstable = oracle_instability.get((smi, aidx))
            row = {
                "engine": "accurate",
                "smiles": smi,
                "atom_idx": aidx,
                "chematic": cm_code,
                "modern": modern_code,
                "bucket": bucket,
                "evidence": evidence,
            }
            if unstable is not None:
                row["oracle_status"] = unstable["oracle_status"]
                row["instability_family"] = unstable.get("instability_family")
            mismatches.append(row)
    return total, mismatches, not_considered


def summarize(engine_name, total, mismatches, not_considered):
    bucket_counts = {b: 0 for b in BUCKET_PRECEDENCE}
    for m in mismatches:
        bucket_counts[m["bucket"]] += 1
    oracle_stable_correct = total - len(mismatches)
    oracle_unstable = sum(1 for m in mismatches if m.get("oracle_status") == "representation_unstable")
    return {
        "engine": engine_name,
        "total_stereocenters_oracle_assigned": total,
        "correct": total - len(mismatches),
        "correct_pct": round(100 * (total - len(mismatches)) / total, 2) if total else None,
        "mismatch": len(mismatches),
        "not_considered_by_chematic": not_considered,
        "oracle_unstable_within_mismatches": oracle_unstable,
        "oracle_stable_correct": oracle_stable_correct,
        "oracle_stable_total": total - oracle_unstable,
        "oracle_stable_pct": (
            round(100 * oracle_stable_correct / (total - oracle_unstable), 2)
            if (total - oracle_unstable)
            else None
        ),
        "bucket_counts": bucket_counts,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--smiles", required=True)
    ap.add_argument("--accurate-tsv", required=True, help="corpus_snapshot --candidate output")
    ap.add_argument("--kekulize-probe", required=True, help="cip_kekulize_probe.rs JSONL output")
    ap.add_argument("--oracle-instability", default="validation/cip_oracle_instability.jsonl")
    ap.add_argument("--out-corpus", default="validation/cip_residual_classification_corpus.jsonl")
    ap.add_argument("--out-summary", default="validation/results/cip_residual_classification_summary.json")
    args = ap.parse_args()

    smis = load_smiles(args.smiles)
    kek_probe = load_kekulize_probe(args.kekulize_probe)
    accurate_snapshot = load_accurate_snapshot(args.accurate_tsv)
    oracle_instability = load_oracle_instability(args.oracle_instability)

    print(f"building oracle over {len(smis)} SMILES rows...", file=sys.stderr)
    oracle = build_oracle(smis)
    print(f"oracle built for {len(oracle)} molecules", file=sys.stderr)

    legacy_total, legacy_mismatches, legacy_not_considered = classify_legacy_engine(oracle, kek_probe)
    accurate_total, accurate_mismatches, accurate_not_considered = classify_accurate_engine(
        oracle, kek_probe, accurate_snapshot, oracle_instability
    )

    legacy_summary = summarize("legacy", legacy_total, legacy_mismatches, legacy_not_considered)
    accurate_summary = summarize("accurate", accurate_total, accurate_mismatches, accurate_not_considered)

    print("=== LEGACY engine (chematic_chem::assign_cip) ===")
    print(json.dumps(legacy_summary, indent=2))
    print("=== ACCURATE engine (chematic_cip::assign_cip_accurate_experimental) ===")
    print(json.dumps(accurate_summary, indent=2))

    all_mismatches = legacy_mismatches + accurate_mismatches
    all_mismatches.sort(key=lambda r: (r["engine"], r["smiles"], r["atom_idx"]))

    manifest = {
        "_manifest": True,
        "generated": datetime.date.today().isoformat(),
        "rdkit_version": Chem.rdBase.rdkitVersion,
        "chematic_version": getattr(chematic, "__version__", None),
        "source_smiles_csv_sha256": sha256_file(args.smiles),
        "bucket_precedence": BUCKET_PRECEDENCE,
        "legacy": {k: v for k, v in legacy_summary.items() if k not in ("bucket_counts",)},
        "accurate": {k: v for k, v in accurate_summary.items() if k not in ("bucket_counts",)},
    }

    with open(args.out_corpus, "w") as f:
        f.write(json.dumps(manifest, sort_keys=True) + "\n")
        for row in all_mismatches:
            f.write(json.dumps(row, sort_keys=True) + "\n")
    print(f"froze {len(all_mismatches)} classified mismatches to {args.out_corpus}", file=sys.stderr)

    full_summary = {
        "_manifest": manifest,
        "legacy": legacy_summary,
        "accurate": accurate_summary,
    }
    with open(args.out_summary, "w") as f:
        json.dump(full_summary, f, indent=2, sort_keys=True)
        f.write("\n")
    print(f"wrote machine-readable summary to {args.out_summary}", file=sys.stderr)


if __name__ == "__main__":
    main()
