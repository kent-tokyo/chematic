#!/usr/bin/env python3
"""Wave 2C (issue #149) audit: ring-constrained E/Z residual diagnosis.

Follows on from PR #229 (Wave 2B), which split the 18 `EZ_SHARED_CANDIDATE_
BOND_RESIDUALS` fixtures into 10 `EZ_SHARED_CARRIER_FULLY_RESOLVED` and 8
`EZ_SHARED_CARRIER_RING_CONSTRAINED_RESIDUALS` (see
`crates/chematic-smiles/src/canonical.rs`, right above the latter constant).
The doc comment there proposes a hypothesis: every one of the 8 residuals'
coupled components includes an alkene end whose own C=X double bond is
*endocyclic* in a 5- or 6-membered ring -- real-world geometry fixed by the
ring, not a free stereochemical choice -- which `compute_stereo_alkene_ends`
has no gate for.

This script tests that hypothesis empirically, matching
`scripts/canonical_residual_diagnosis.py`'s conventions (buckets measured
SEPARATELY not pooled, self-test, JSONL + summary-JSON output, live RDKit
oracle -- nothing assumed). Diagnosis only: no production code is touched by
this script or by the Rust example it drives
(`crates/chematic-smiles/examples/ez_ring_constrained_residual_audit.rs`).

Two independent chematic-side data pulls (both via `cargo run`, never
`import chematic` -- this avoids depending on whatever the active venv's
installed bindings happen to be built from, matching this repo's
oracle-must-be-independent convention):
  1. The 18 pinned fixtures (no corpus arg) -- small enough to classify each
     one individually, never pooled (per this repo's own established norm:
     "never silently pool the 8 residual fixtures' results into an aggregate
     that could hide a fixture where the hypothesis doesn't actually hold").
  2. The full 5,000-molecule committed corpus
     (`scripts/descriptor_census_corpus.smi`) -- for the corpus-wide
     blast-radius measurement, which needs the full stereo-alkene-end
     population, not just the 31 atoms in a size-2 coupling component.

For every row (one stereo-alkene end), independently cross-checked against
live RDKit (`Chem.FindPotentialStereo` + `Chem.AssignStereochemistry`) with
an EXPLICIT atom-index-correspondence check first (chematic and RDKit both
assign heavy-atom indices in SMILES-encounter order for the same input
string, but this is verified per row via element-symbol match, never
assumed -- rows that fail are flagged and excluded from RDKit-oracle-based
counts, not silently included with a guessed verdict).

Three candidate gating-rule blast radii are measured INDIVIDUALLY, never
combined into one number (per this repo's established norm: ring-size
threshold, RDKit-not-potential-stereo, and ring-topology-specific
impossibility are three separate empirical questions):
  (a) end atom's own ring size < N -- measured BOTH as "the end atom sits in
      a small ring" (over-inclusive: also flags exocyclic-but-real ends that
      merely happen to sit in a small ring, like atom 16 in the fixture-1
      worked example) and as "the end atom's OWN double bond is endocyclic
      in a small ring" (the actual predicate the hypothesis needs).
  (b) RDKit does not list this bond in `FindPotentialStereo`'s `Bond_Double`
      output at all (`cleanIt=False, flagPossible=True`) -- this RDKit
      version (2026.03.3) has no `NOT_POSSIBLE` value in `StereoSpecified`
      (only Unspecified/Specified/Unknown -- confirmed by direct enum dump,
      see `docs/ez_ring_constrained_residual_audit.md`), so "absent from the
      list" is the operational equivalent, confirmed mechanistically:
      1,2-disubstituted cyclohexene is absent, 1,2-disubstituted cyclooctene
      (large enough for real trans) is present as Unspecified.
  (c) same absence, but RESTRICTED to bonds chematic itself found endocyclic
      -- genuinely narrower than (b) only if (b) also catches non-ring
      exclusions (e.g. symmetric-substituent cases); if (b) and (c) produce
      the SAME count on this corpus, that is reported as a finding (a
      coincidence, not manufactured as a third distinct number).

"ends excluded" under a rule is directly measurable without touching
production code; "canonical output would actually change" is NOT (that
requires implementing the rule). This script reports excluded-and-currently-
marked ends as an explicit UPPER BOUND on output change, never as a
measurement of output change itself.

Usage:
    .venv/bin/python3 scripts/ez_ring_constrained_residual_diagnosis.py [CORPUS]
    .venv/bin/python3 scripts/ez_ring_constrained_residual_diagnosis.py --self-test

Writes:
    validation/results/ez_ring_constrained_residual_audit.jsonl          (per-end rows, full corpus)
    validation/results/ez_ring_constrained_residual_audit_summary.json   (fixture classification + blast radius)
"""
import argparse
import json
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_JSONL = os.path.join(ROOT, "validation", "results", "ez_ring_constrained_residual_audit.jsonl")
OUT_SUMMARY = os.path.join(
    ROOT, "validation", "results", "ez_ring_constrained_residual_audit_summary.json"
)
DEFAULT_CORPUS = os.path.join(ROOT, "scripts", "descriptor_census_corpus.smi")

# Mirrors `canonical.rs`'s two constants exactly -- used ONLY to label which
# of the 18 pinned fixtures each row came from; never re-derived, never
# hand-edited independently of the source of truth in canonical.rs.
EZ_SHARED_CARRIER_FULLY_RESOLVED = [
    r"CCCCC/N=c1\c(O)c(O)\c1=N/[C@@H](Cc1ccc(NC(=O)c2c(Cl)cncc2Cl)cc1)C(=O)O",
    r"O=C(Nc1ccc(C[C@H](/N=c2\c(O)c(O)\c2=N/Cc2ccccc2)C(=O)O)cc1)c1c(Cl)cncc1Cl",
    r"CCC/N=c1\c(O)c(O)\c1=N/[C@@H](Cc1ccc(NC(=O)c2c(Cl)cncc2Cl)cc1)C(=O)O",
    r"O=C(Nc1ccc(C[C@H](/N=c2\c(O)c(O)\c2=N/c2ccccc2)C(=O)O)cc1)c1c(Cl)cncc1Cl",
    r"CC(C)(C)/N=c1\c(O)c(O)\c1=N/[C@@H](Cc1ccc(NC(=O)c2c(Cl)cncc2Cl)cc1)C(=O)O",
    r"CCOC(=O)C1=C(C)N=C(C)/C(=C(/O)OCC)C1c1cccc(I)c1",
    r"CCO/C(O)=C1\C(C)=NC(C)=C(C(=O)OC)C1c1ccccc1C(F)(F)F",
    r"CCOC(=O)C1=C(C)N=C(C)/C(=C(/O)OCC)C1c1ccccc1OC",
    r"CCO/C(O)=C1\C(C)=NC(C)=C(C(=O)OC)C1c1ccc([N+](=O)[O-])cc1",
    r"CCO/C(O)=C1\C(C)=NC(C)=C(C(=O)OC)C1c1cccc(C(F)(F)F)c1",
]

EZ_SHARED_CARRIER_RING_CONSTRAINED_RESIDUALS = [
    r"CC1=C2CC[C@H](/C=N/N=C(N)N)[C@@]2(C)CC/C1=N\N=C(N)N",
    r"CC1=C2CC[C@@H](/C=N/N=C(N)N)[C@@]2(C)CC/C1=N\N=C(N)N",
    r"COC(=O)/C=C/[C@H]1CCC2=C(C)/C(=N/N=C(N)N)CC[C@@]21C",
    r"CCOC(=O)C1=C(C)N=C(C)/C(=C(/O)OCC)C1c1ccc(I)cc1",
    r"CCO/C(O)=C1\C(C)=NC(C)=C(C(=O)OC)C1c1ccccc1[N+](=O)[O-]",
    r"CCOC(=O)C1=C(C)N=C(C)/C(=C(/O)OCC)C1c1cccc(C)c1",
    r"CCOC(=O)C1=C(C)N=C(C)/C(=C(/O)OCC)C1c1cccc(OC)c1",
    r"CCO/C(O)=C(\C1=NCCN1)c1nnc(N)s1",
]

RING_SIZE_THRESHOLDS = (7, 8)  # excludes rings of size < threshold


def fixture_label(smiles):
    if smiles in EZ_SHARED_CARRIER_FULLY_RESOLVED:
        return "fully_resolved"
    if smiles in EZ_SHARED_CARRIER_RING_CONSTRAINED_RESIDUALS:
        return "ring_constrained_residual"
    return "corpus"


# ---------------------------------------------------------------------------
# Chematic side: drive the Rust example, never `import chematic` (see module
# docstring -- avoids depending on whichever bindings happen to be installed).
# ---------------------------------------------------------------------------

def run_chematic_example(corpus_path=None):
    cmd = [
        "cargo", "run", "-p", "chematic-smiles", "--release",
        "--example", "ez_ring_constrained_residual_audit",
    ]
    if corpus_path is not None:
        cmd += ["--", corpus_path]
    result = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"cargo run failed (exit {result.returncode}):\n{result.stderr[-4000:]}"
        )
    rows = []
    for line in result.stdout.splitlines():
        line = line.strip()
        if not line:
            continue
        rows.append(json.loads(line))
    return rows


# ---------------------------------------------------------------------------
# RDKit side: independent per-molecule oracle. Explicit index-correspondence
# check FIRST (element-symbol match on both the end atom and its double-bond
# partner) -- never assumed just because both parsers happen to usually
# agree on SMILES-encounter-order indexing.
# ---------------------------------------------------------------------------

def rdkit_molecule_oracle(smiles):
    """Returns (rdkit_mol_or_None, potential_by_bond_idx, assigned_by_bond_idx).

    potential_by_bond_idx: {bond_idx: (specified_str, descriptor_str)} for
    every Bond_Double entry `Chem.FindPotentialStereo` reports (absence from
    this dict means RDKit did not even consider the bond potentially
    stereogenic).
    assigned_by_bond_idx: {bond_idx: bond_stereo_str} from
    `Chem.AssignStereochemistry(cleanIt=True, force=True)`, i.e. RDKit's
    real E/Z semantic label (or STEREONONE/STEREOANY) for the ORIGINAL input
    geometry.
    """
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, {}, {}

    potential = {}
    for si in Chem.FindPotentialStereo(mol, cleanIt=False, flagPossible=True):
        if str(si.type) == "Bond_Double":
            potential[si.centeredOn] = (str(si.specified), str(si.descriptor))

    assigned_mol = Chem.Mol(mol)
    Chem.AssignStereochemistry(assigned_mol, cleanIt=True, force=True)
    assigned = {}
    for b in assigned_mol.GetBonds():
        if b.GetBondType() == Chem.BondType.DOUBLE:
            assigned[b.GetIdx()] = str(b.GetStereo())

    return mol, potential, assigned


def crosscheck_row(rdkit_mol, potential, assigned, row):
    """Augments one chematic-side row with RDKit-oracle fields, verifying
    atom-index correspondence explicitly before trusting any RDKit verdict."""
    out = dict(row)
    end_idx = row["end_atom_idx"]
    partner_idx = row["partner_atom_idx"]

    if rdkit_mol is None:
        out["rdkit_index_correspondence_ok"] = False
        out["rdkit_correspondence_failure_reason"] = "rdkit_parse_failed"
        return _fill_rdkit_nulls(out)

    n = rdkit_mol.GetNumAtoms()
    if end_idx >= n or partner_idx >= n:
        out["rdkit_index_correspondence_ok"] = False
        out["rdkit_correspondence_failure_reason"] = "index_out_of_range"
        return _fill_rdkit_nulls(out)

    end_atom = rdkit_mol.GetAtomWithIdx(end_idx)
    partner_atom = rdkit_mol.GetAtomWithIdx(partner_idx)
    if end_atom.GetSymbol() != row["end_element"] or partner_atom.GetSymbol() != row["partner_element"]:
        out["rdkit_index_correspondence_ok"] = False
        out["rdkit_correspondence_failure_reason"] = "element_mismatch"
        return _fill_rdkit_nulls(out)

    bond = rdkit_mol.GetBondBetweenAtoms(end_idx, partner_idx)
    from rdkit import Chem
    if bond is None or bond.GetBondType() != Chem.BondType.DOUBLE:
        out["rdkit_index_correspondence_ok"] = False
        out["rdkit_correspondence_failure_reason"] = "no_matching_double_bond"
        return _fill_rdkit_nulls(out)

    out["rdkit_index_correspondence_ok"] = True
    out["rdkit_correspondence_failure_reason"] = None
    bond_idx = bond.GetIdx()
    out["rdkit_bond_idx"] = bond_idx

    pot = potential.get(bond_idx)
    out["rdkit_potential_stereo"] = pot is not None
    out["rdkit_specified"] = pot[0] if pot else None
    out["rdkit_descriptor"] = pot[1] if pot else None

    stereo = assigned.get(bond_idx, "STEREONONE")
    out["rdkit_bond_stereo"] = stereo
    out["rdkit_ez_assignable"] = stereo not in ("STEREONONE", "STEREOANY")
    return out


def _fill_rdkit_nulls(out):
    out["rdkit_bond_idx"] = None
    out["rdkit_potential_stereo"] = None
    out["rdkit_specified"] = None
    out["rdkit_descriptor"] = None
    out["rdkit_bond_stereo"] = None
    out["rdkit_ez_assignable"] = None
    return out


def crosscheck_rows(rows):
    """Groups rows by SMILES so the RDKit oracle runs once per molecule, not
    once per end -- matching `canonical_residual_diagnosis.py`'s per-molecule
    granularity."""
    by_smiles = {}
    for row in rows:
        by_smiles.setdefault(row["smiles"], []).append(row)

    out_rows = []
    for smiles, group in by_smiles.items():
        rdkit_mol, potential, assigned = rdkit_molecule_oracle(smiles)
        for row in group:
            out_rows.append(crosscheck_row(rdkit_mol, potential, assigned, row))
    return out_rows


# ---------------------------------------------------------------------------
# Blast-radius rules -- each evaluated INDEPENDENTLY (never combined).
# ---------------------------------------------------------------------------

def rule_a_atom_in_small_ring(row, threshold):
    return bool(row["end_atom_ring_sizes"]) and min(row["end_atom_ring_sizes"]) < threshold


def rule_a_bond_endocyclic_small_ring(row, threshold):
    return row["double_bond_endocyclic"] and min(row["double_bond_endocyclic_ring_sizes"]) < threshold


def rule_b_rdkit_not_potential(row):
    if not row["rdkit_index_correspondence_ok"]:
        return None  # unmeasurable for this row
    return row["rdkit_potential_stereo"] is False


def rule_c_ring_topology_impossible(row):
    b = rule_b_rdkit_not_potential(row)
    if b is None:
        return None
    return b and row["double_bond_endocyclic"]


RULES = [
    ("a_atom_in_ring_lt7", lambda r: rule_a_atom_in_small_ring(r, 7)),
    ("a_atom_in_ring_lt8", lambda r: rule_a_atom_in_small_ring(r, 8)),
    ("a_bond_endocyclic_lt7", lambda r: rule_a_bond_endocyclic_small_ring(r, 7)),
    ("a_bond_endocyclic_lt8", lambda r: rule_a_bond_endocyclic_small_ring(r, 8)),
    ("b_rdkit_not_potential_stereo", rule_b_rdkit_not_potential),
    ("c_ring_topology_impossible", rule_c_ring_topology_impossible),
]


def blast_radius_table(rows):
    """Per rule: ends excluded (count/fraction, unmeasurable rows tracked
    separately, never silently dropped from the denominator), how many of
    those currently carry a marker (upper bound on output change -- NOT a
    measurement of output change), how many are coupled, and distinct
    molecule counts."""
    total_ends = len(rows)
    table = {}
    for name, rule_fn in RULES:
        excluded_rows = []
        unmeasurable = 0
        for row in rows:
            verdict = rule_fn(row)
            if verdict is None:
                unmeasurable += 1
                continue
            if verdict:
                excluded_rows.append(row)

        excluded_marked = [
            r for r in excluded_rows
            if r.get("marker_placed") is True
        ]
        excluded_coupled = [r for r in excluded_rows if r["coupled"]]
        molecules_excluded = {r["smiles"] for r in excluded_rows}
        molecules_excluded_and_marked = {r["smiles"] for r in excluded_marked}

        table[name] = {
            "total_ends": total_ends,
            "unmeasurable_ends": unmeasurable,
            "ends_excluded": len(excluded_rows),
            "ends_excluded_pct_of_measurable": round(
                100 * len(excluded_rows) / (total_ends - unmeasurable), 2
            ) if (total_ends - unmeasurable) else None,
            "excluded_and_currently_marked_UPPER_BOUND_on_output_change": len(excluded_marked),
            "excluded_and_coupled": len(excluded_coupled),
            "distinct_molecules_with_excluded_end": len(molecules_excluded),
            "distinct_molecules_with_excluded_and_marked_end": len(molecules_excluded_and_marked),
        }
    return table


def per_fixture_rule_effect(fixture_rows_by_smiles, fixture_list, rule_name, rule_fn):
    """Per pinned fixture: does this rule exclude >=1 end, and >=1
    excluded-and-marked end -- reported per fixture, not pooled, per this
    repo's own established norm."""
    out = []
    for smiles in fixture_list:
        rows = fixture_rows_by_smiles.get(smiles, [])
        excluded = []
        for row in rows:
            v = rule_fn(row)
            if v:
                excluded.append(row)
        excluded_marked = [r for r in excluded if r.get("marker_placed") is True]
        out.append({
            "smiles": smiles,
            "n_ends": len(rows),
            "n_excluded": len(excluded),
            "n_excluded_and_marked": len(excluded_marked),
            "excluded_end_atom_idxs": sorted(r["end_atom_idx"] for r in excluded),
        })
    return out


# ---------------------------------------------------------------------------
# Per-fixture individual classification (never pooled) for the 18 pinned
# fixtures.
# ---------------------------------------------------------------------------

def classify_fixtures(fixture_rows):
    by_smiles = {}
    for row in fixture_rows:
        by_smiles.setdefault(row["smiles"], []).append(row)

    classification = []
    for smiles in EZ_SHARED_CARRIER_FULLY_RESOLVED + EZ_SHARED_CARRIER_RING_CONSTRAINED_RESIDUALS:
        rows = by_smiles.get(smiles, [])
        coupled_rows = [r for r in rows if r["coupled"]]
        classification.append({
            "smiles": smiles,
            "label": fixture_label(smiles),
            "n_stereo_alkene_ends": len(rows),
            "n_coupled_ends": len(coupled_rows),
            "coupled_ends": [
                {
                    "end_atom_idx": r["end_atom_idx"],
                    "end_element": r["end_element"],
                    "partner_atom_idx": r["partner_atom_idx"],
                    "partner_element": r["partner_element"],
                    "component_other_members": r["component_other_members"],
                    "double_bond_endocyclic": r["double_bond_endocyclic"],
                    "double_bond_endocyclic_ring_sizes": r["double_bond_endocyclic_ring_sizes"],
                    "end_atom_ring_sizes": r["end_atom_ring_sizes"],
                    "marker_placed": r.get("marker_placed"),
                    "rdkit_index_correspondence_ok": r.get("rdkit_index_correspondence_ok"),
                    "rdkit_potential_stereo": r.get("rdkit_potential_stereo"),
                    "rdkit_specified": r.get("rdkit_specified"),
                    "rdkit_bond_stereo": r.get("rdkit_bond_stereo"),
                }
                for r in coupled_rows
            ],
            # The hypothesis's own predicate: at least one end of the
            # coupled pair has an endocyclic double bond AND RDKit
            # independently confirms that specific bond is not a real
            # stereocenter (STEREONONE after AssignStereochemistry).
            "hypothesis_holds": any(
                r["double_bond_endocyclic"] and r.get("rdkit_bond_stereo") == "STEREONONE"
                for r in coupled_rows
            ) if coupled_rows else None,
        })
    return classification


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run(corpus_path, out_jsonl=OUT_JSONL, out_summary=OUT_SUMMARY):
    print(f"running chematic example on 18 pinned fixtures...", file=sys.stderr)
    fixture_rows_raw = run_chematic_example(corpus_path=None)
    fixture_rows = crosscheck_rows(fixture_rows_raw)

    print(f"running chematic example on corpus {corpus_path} ...", file=sys.stderr)
    corpus_rows_raw = run_chematic_example(corpus_path=corpus_path)
    corpus_rows = crosscheck_rows(corpus_rows_raw)

    os.makedirs(os.path.dirname(out_jsonl), exist_ok=True)
    with open(out_jsonl, "w") as f:
        for row in corpus_rows:
            f.write(json.dumps(row) + "\n")

    classification = classify_fixtures(fixture_rows)

    blast_radius = blast_radius_table(corpus_rows)

    fixture_rows_by_smiles = {}
    for row in fixture_rows:
        fixture_rows_by_smiles.setdefault(row["smiles"], []).append(row)
    per_fixture_effects = {}
    for name, rule_fn in RULES:
        per_fixture_effects[name] = {
            "fully_resolved": per_fixture_rule_effect(
                fixture_rows_by_smiles, EZ_SHARED_CARRIER_FULLY_RESOLVED, name, rule_fn
            ),
            "ring_constrained_residual": per_fixture_rule_effect(
                fixture_rows_by_smiles, EZ_SHARED_CARRIER_RING_CONSTRAINED_RESIDUALS, name, rule_fn
            ),
        }

    n_correspondence_fail_corpus = sum(
        1 for r in corpus_rows if not r["correspondence_ok"]
    )
    n_rdkit_correspondence_fail_corpus = sum(
        1 for r in corpus_rows if not r["rdkit_index_correspondence_ok"]
    )

    summary = {
        "corpus": corpus_path,
        "n_fixture_rows": len(fixture_rows),
        "n_corpus_rows": len(corpus_rows),
        "n_corpus_molecules_with_ends": len({r["smiles"] for r in corpus_rows}),
        "n_corpus_coupled_ends": sum(1 for r in corpus_rows if r["coupled"]),
        "chematic_correspondence_failures_corpus": n_correspondence_fail_corpus,
        "rdkit_index_correspondence_failures_corpus": n_rdkit_correspondence_fail_corpus,
        "fixture_classification": classification,
        "blast_radius": blast_radius,
        "blast_radius_per_fixture_effect": per_fixture_effects,
    }
    with open(out_summary, "w") as f:
        json.dump(summary, f, indent=2)

    print(json.dumps({k: v for k, v in summary.items() if k not in (
        "fixture_classification", "blast_radius_per_fixture_effect"
    )}, indent=2))
    print(f"\nwrote {len(corpus_rows)} corpus rows -> {os.path.relpath(out_jsonl, ROOT)}")
    print(f"wrote summary -> {os.path.relpath(out_summary, ROOT)}")
    return summary


# ---------------------------------------------------------------------------
# Self-test: verifies the rule/classification functions actually
# discriminate (fail-closed, positive AND negative controls -- per this
# repo's measurement-harness-controls convention).
# ---------------------------------------------------------------------------

def self_test():
    row_base = {
        "smiles": "X", "end_atom_idx": 0, "end_element": "C",
        "partner_atom_idx": 1, "partner_element": "N",
        "coupled": True, "component_other_members": [2],
        "end_atom_ring_sizes": [6], "double_bond_endocyclic": True,
        "double_bond_endocyclic_ring_sizes": [6],
        "correspondence_ok": True, "marker_placed": False,
        "rdkit_index_correspondence_ok": True,
        "rdkit_potential_stereo": False,
    }

    def with_(**kw):
        d = dict(row_base)
        d.update(kw)
        return d

    # rule a: atom-in-ring vs bond-endocyclic must differ on the "atom in
    # small ring but OWN double bond exocyclic" case (fixture-1 atom 16
    # shape).
    exocyclic_but_ring_atom = with_(
        end_atom_ring_sizes=[6], double_bond_endocyclic=False,
        double_bond_endocyclic_ring_sizes=[],
    )
    assert rule_a_atom_in_small_ring(exocyclic_but_ring_atom, 7) is True
    assert rule_a_bond_endocyclic_small_ring(exocyclic_but_ring_atom, 7) is False
    print("rule a atom-vs-bond distinction: OK (exocyclic-in-small-ring case differs)")

    endocyclic_6ring = with_()
    assert rule_a_bond_endocyclic_small_ring(endocyclic_6ring, 7) is True
    assert rule_a_bond_endocyclic_small_ring(endocyclic_6ring, 6) is False
    print("rule a threshold sensitivity: OK")

    endocyclic_8ring = with_(
        double_bond_endocyclic=True, double_bond_endocyclic_ring_sizes=[8],
        end_atom_ring_sizes=[8],
    )
    assert rule_a_bond_endocyclic_small_ring(endocyclic_8ring, 7) is False
    assert rule_a_bond_endocyclic_small_ring(endocyclic_8ring, 8) is False
    print("rule a large-ring negative control: OK (8-ring not excluded at threshold 7 or 8)")

    # rule b/c
    not_potential = with_(rdkit_index_correspondence_ok=True, rdkit_potential_stereo=False)
    assert rule_b_rdkit_not_potential(not_potential) is True
    assert rule_c_ring_topology_impossible(not_potential) is True  # also endocyclic in row_base

    potential = with_(rdkit_potential_stereo=True)
    assert rule_b_rdkit_not_potential(potential) is False
    assert rule_c_ring_topology_impossible(potential) is False

    exocyclic_not_potential = with_(
        rdkit_potential_stereo=False, double_bond_endocyclic=False,
        double_bond_endocyclic_ring_sizes=[],
    )
    assert rule_b_rdkit_not_potential(exocyclic_not_potential) is True
    assert rule_c_ring_topology_impossible(exocyclic_not_potential) is False, (
        "rule c must be narrower than rule b: an exocyclic bond RDKit "
        "excludes for a non-ring reason must not count as ring-topology-caused"
    )
    print("rule b/c distinction: OK (exocyclic RDKit-exclusion does not trip rule c)")

    unmeasurable = with_(rdkit_index_correspondence_ok=False)
    assert rule_b_rdkit_not_potential(unmeasurable) is None
    assert rule_c_ring_topology_impossible(unmeasurable) is None
    print("unmeasurable-row handling: OK (correspondence failure -> None, not guessed)")

    # blast_radius_table denominator/unmeasurable accounting
    rows = [row_base, potential, unmeasurable]
    table = blast_radius_table(rows)
    assert table["b_rdkit_not_potential_stereo"]["total_ends"] == 3
    assert table["b_rdkit_not_potential_stereo"]["unmeasurable_ends"] == 1
    assert table["b_rdkit_not_potential_stereo"]["ends_excluded"] == 1  # row_base only
    print("blast_radius_table accounting: OK")

    # crosscheck_row correspondence-failure paths, using a live RDKit mol
    from rdkit import Chem
    mol = Chem.MolFromSmiles("C=C")
    potential_map, assigned_map = {}, {}
    bad_element_row = {
        "smiles": "C=C", "end_atom_idx": 0, "end_element": "N",  # wrong on purpose
        "partner_atom_idx": 1, "partner_element": "C",
    }
    out = crosscheck_row(mol, potential_map, assigned_map, bad_element_row)
    assert out["rdkit_index_correspondence_ok"] is False
    assert out["rdkit_correspondence_failure_reason"] == "element_mismatch"
    assert out["rdkit_potential_stereo"] is None
    print("crosscheck_row element-mismatch detection: OK")

    good_row = {
        "smiles": "C=C", "end_atom_idx": 0, "end_element": "C",
        "partner_atom_idx": 1, "partner_element": "C",
    }
    out = crosscheck_row(mol, potential_map, assigned_map, good_row)
    assert out["rdkit_index_correspondence_ok"] is True
    assert out["rdkit_potential_stereo"] is False  # ethylene has no real E/Z (H substituents only)
    print("crosscheck_row good-correspondence path: OK")

    # Positive control: fixture 1's known culprit (atom 1, ring-endocyclic
    # C1=C2, coupled with atom 16) must show rdkit_bond_stereo ==
    # STEREONONE and rdkit_potential_stereo == False -- this is the exact
    # RDKit-side confirmation the classification/verdict depends on. If this
    # ever stops reproducing, the whole audit's empirical basis needs
    # re-checking, not just this test.
    smi = r"CC1=C2CC[C@H](/C=N/N=C(N)N)[C@@]2(C)CC/C1=N\N=C(N)N"
    rdkit_mol, potential, assigned = rdkit_molecule_oracle(smi)
    row1 = {"smiles": smi, "end_atom_idx": 1, "end_element": "C",
            "partner_atom_idx": 2, "partner_element": "C"}
    out1 = crosscheck_row(rdkit_mol, potential, assigned, row1)
    assert out1["rdkit_index_correspondence_ok"] is True
    assert out1["rdkit_bond_stereo"] == "STEREONONE", out1
    assert out1["rdkit_potential_stereo"] is False, out1
    row16 = {"smiles": smi, "end_atom_idx": 16, "end_element": "C",
             "partner_atom_idx": 17, "partner_element": "N"}
    out16 = crosscheck_row(rdkit_mol, potential, assigned, row16)
    assert out16["rdkit_bond_stereo"] == "STEREOE", out16
    assert out16["rdkit_potential_stereo"] is True, out16
    print("positive control (fixture-1 atom1/atom16 RDKit split): OK -- "
          "reproduces the culprit/free-partner asymmetry the audit depends on")

    print("\nself-test OK")
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("corpus", nargs="?", default=DEFAULT_CORPUS)
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--out-jsonl", default=OUT_JSONL)
    ap.add_argument("--out-summary", default=OUT_SUMMARY)
    args = ap.parse_args()

    if args.self_test:
        return self_test()

    run(args.corpus, args.out_jsonl, args.out_summary)
    return 0


if __name__ == "__main__":
    sys.exit(main())
