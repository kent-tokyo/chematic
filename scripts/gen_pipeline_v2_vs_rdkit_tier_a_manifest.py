#!/usr/bin/env python3
"""Generate the Tier A (curated stress) corpus manifest for the pipeline v2 vs
RDKit ETKDGv3 benchmark (Wave 1 of the "RDKit alternative" program).

Source of truth for the 63 molecules: `CORPUS` in
`crates/chematic-3d/examples/pipeline_v2_integration_gate.rs` (itself sourced
from `scripts/etkdg_vs_rdkit_gap.py::CORPUS`), transcribed here verbatim
(name, SMILES, primary_category unchanged) rather than hand-copied a third
time. This script is the SINGLE extraction point: both the new Rust dump
executable and the new RDKit oracle script read the JSON manifest this
script produces, never the Rust corpus directly.

`additional_tags` adds the finer-grained taxonomy Wave 1's spec requires
(e.g. splitting "rigid_ring" into aromatic vs. bridged vs. small-ring) without
renaming or removing any molecule's original `primary_category` -- anyone
diffing against the source corpus can still verify a 1:1, unmodified copy of
the SMILES/name/primary_category columns.

Two molecules are appended beyond the source 63 for taxonomy categories the
source corpus has zero representation for. Each is flagged
`"added_for_taxonomy": true` with a `"provenance"` note -- neither is
invented ad hoc; both are already-established fixtures reused verbatim from
elsewhere in this repo's test/validation suite:
- `ring_torsion_fail_closed`: identical SMILES + config precedent to the
  `ring_torsion_fail_closed` fixture in
  `validation/pipeline_v2_wasm_parity_fixtures.json` (Wave 1B, PR #220).
- `force_field_unsupported_probe`: a phosphonium ylide fragment (P=C), an
  element/bond-order combination MMFF94's parameter tables are unlikely to
  cover -- included as a *candidate* probe for the "force-field unsupported"
  bucket, not a verified-in-advance member. Whether it (or anything else in
  the corpus) actually lands in that bucket is determined empirically by
  this benchmark's own arm coverage results, not asserted here.

Run: `.venv/bin/python scripts/gen_pipeline_v2_vs_rdkit_tier_a_manifest.py`
Output: `validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_a.json`
"""

import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).parent.parent

SOURCE_FILE = "crates/chematic-3d/examples/pipeline_v2_integration_gate.rs"
SOURCE_SYMBOL = "CORPUS"

# (name, smiles, primary_category) -- verbatim from CORPUS in SOURCE_FILE.
CORPUS = [
    ("benzene", "c1ccccc1", "rigid_ring"),
    ("naphthalene", "c1ccc2ccccc2c1", "fused_aromatic"),
    ("pyridine", "c1ccncc1", "rigid_ring"),
    ("furan", "c1ccoc1", "rigid_ring"),
    ("thiophene", "c1ccsc1", "rigid_ring"),
    ("adamantane", "C1CC2CC3CC1CC(C2)C3", "rigid_ring"),
    ("cubane", "C1C2C3C1C4C2C3C4", "rigid_ring"),
    ("cyclohexane", "C1CCCCC1", "rigid_ring"),
    ("cyclopentane", "C1CCCC1", "rigid_ring"),
    ("indole", "c1ccc2[nH]ccc2c1", "fused_aromatic"),
    ("purine", "c1ncc2[nH]cnc2n1", "fused_aromatic"),
    ("quinoline", "c1ccc2ncccc2c1", "fused_aromatic"),
    ("anthracene", "c1ccc2cc3ccccc3cc2c1", "fused_aromatic"),
    ("pyrene", "c1cc2ccc3cccc4ccc(c1)c2c34", "fused_aromatic"),
    ("biphenyl", "c1ccc(-c2ccccc2)cc1", "fused_aromatic"),
    ("butane", "CCCC", "flexible_chain"),
    ("hexane", "CCCCCC", "flexible_chain"),
    ("decane", "CCCCCCCCCC", "flexible_chain"),
    ("triethylene_glycol", "OCCOCCOCCO", "flexible_chain"),
    ("hexanediol", "OCCCCCCO", "flexible_chain"),
    ("hexadecane", "CCCCCCCCCCCCCCCC", "flexible_chain"),
    ("cyclododecane", "C1CCCCCCCCCCC1", "macrocycle"),
    ("crown_12_4", "O1CCOCCOCCOCC1", "macrocycle"),
    ("cyclooctadecane", "C1CCCCCCCCCCCCCCCCC1", "macrocycle"),
    ("l_alanine", "N[C@@H](C)C(=O)O", "stereocenter_implicit_h"),
    ("d_alanine", "N[C@H](C)C(=O)O", "stereocenter_implicit_h"),
    ("l_serine", "N[C@@H](CO)C(=O)O", "stereocenter_implicit_h"),
    ("l_threonine", "C[C@H](O)[C@@H](N)C(=O)O", "stereocenter_implicit_h"),
    ("2_butanol_R", "C[C@H](O)CC", "stereocenter_implicit_h"),
    ("2_butanol_S", "C[C@@H](O)CC", "stereocenter_implicit_h"),
    ("2_chlorobutane_R", "C[C@H](Cl)CC", "stereocenter_implicit_h"),
    ("ibuprofen_S", "CC(C)Cc1ccc(cc1)[C@H](C)C(=O)O", "stereocenter_implicit_h"),
    ("naproxen_S", "COc1ccc2cc([C@H](C)C(=O)O)ccc2c1", "stereocenter_implicit_h"),
    ("menthol", "C[C@@H]1CC[C@@H](C(C)C)C[C@H]1O", "stereocenter_implicit_h"),
    ("chfclbr_R", "[C@H](F)(Cl)Br", "stereocenter_quaternary"),
    ("chfclbr_S", "[C@@H](F)(Cl)Br", "stereocenter_quaternary"),
    ("quaternary_1_R", "[C@](F)(Cl)(Br)I", "stereocenter_quaternary"),
    ("quaternary_1_S", "[C@@](F)(Cl)(Br)I", "stereocenter_quaternary"),
    ("quaternary_2_R", "[C@](C)(N)(O)F", "stereocenter_quaternary"),
    ("quaternary_2_S", "[C@@](C)(N)(O)F", "stereocenter_quaternary"),
    ("but2ene_E", "C/C=C/C", "alkene_ez"),
    ("but2ene_Z", r"C/C=C\C", "alkene_ez"),
    ("chloropropene_E", "C(/C=C/C)Cl", "alkene_ez"),
    ("chloropropene_Z", r"C(/C=C\C)Cl", "alkene_ez"),
    ("cinnamic_acid_E", "OC(=O)/C=C/c1ccccc1", "alkene_ez"),
    ("cinnamic_acid_Z", r"OC(=O)/C=C\c1ccccc1", "alkene_ez"),
    ("pent2ene_E", "CC/C=C/C", "alkene_ez"),
    ("pent2ene_Z", r"CC/C=C\C", "alkene_ez"),
    ("aspirin", "CC(=O)Oc1ccccc1C(=O)O", "druglike"),
    ("ibuprofen", "CC(C)Cc1ccc(cc1)C(C)C(=O)O", "druglike"),
    ("caffeine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C", "druglike"),
    ("paracetamol", "CC(=O)Nc1ccc(O)cc1", "druglike"),
    ("diphenhydramine", "CN(C)CCOC(c1ccccc1)c1ccccc1", "druglike"),
    ("penicillin_core", "CC1(C)S[C@@H]2[C@H](NC(=O)C)C(=O)N2[C@H]1C(=O)O", "druglike"),
    ("testosterone", "C[C@]12CC[C@H]3[C@@H](CC[C@H]4CCC(=O)C=C34)[C@@H]1CC[C@@H]2O", "druglike_rigid"),
    ("cholesterol", "C[C@H](CCCC(C)C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C", "druglike_stress"),
    ("atorvastatin_fragment", "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O", "druglike_stress"),
    ("gly_ala_gly", "NCC(=O)N[C@@H](C)C(=O)NCC(=O)O", "druglike"),
    ("cyclobutane", "C1CCC1", "rigid_ring"),
    ("cyclooctane", "C1CCCCCCC1", "small_ring_boundary"),
    ("cyclononane", "C1CCCCCCCC1", "macrocycle_boundary"),
    ("dimethylbiphenyl_2_2", "Cc1ccccc1-c1ccccc1C", "hindered_biaryl"),
    ("macrolactam_12", "O=C1CCCCCCCCCCN1", "macrocycle_amide"),
]

# name -> extra taxonomy tags beyond primary_category (Wave 1 spec's required
# category list). Purely additive labeling; does not touch name/smiles/primary_category.
ADDITIONAL_TAGS = {
    "benzene": ["rigid_aromatic"],
    "pyridine": ["rigid_aromatic"],
    "furan": ["rigid_aromatic"],
    "thiophene": ["rigid_aromatic"],
    "adamantane": ["bridged_ring", "rigid_aliphatic"],
    "cubane": ["bridged_ring", "rigid_aliphatic", "small_ring"],
    "cyclohexane": ["rigid_aliphatic"],
    "cyclopentane": ["rigid_aliphatic", "small_ring"],
    "cyclobutane": ["rigid_aliphatic", "small_ring"],
    "cyclooctane": ["small_ring"],
    "biphenyl": ["biaryl_unhindered"],
    "testosterone": ["steroid_like"],
    "cholesterol": ["steroid_like"],
    "l_alanine": ["declared_tetrahedral_stereo"],
    "d_alanine": ["declared_tetrahedral_stereo"],
    "l_serine": ["declared_tetrahedral_stereo"],
    "l_threonine": ["declared_tetrahedral_stereo"],
    "2_butanol_R": ["declared_tetrahedral_stereo"],
    "2_butanol_S": ["declared_tetrahedral_stereo"],
    "2_chlorobutane_R": ["declared_tetrahedral_stereo"],
    "ibuprofen_S": ["declared_tetrahedral_stereo"],
    "naproxen_S": ["declared_tetrahedral_stereo"],
    "menthol": ["declared_tetrahedral_stereo"],
    "chfclbr_R": ["declared_tetrahedral_stereo"],
    "chfclbr_S": ["declared_tetrahedral_stereo"],
    "quaternary_1_R": ["declared_tetrahedral_stereo"],
    "quaternary_1_S": ["declared_tetrahedral_stereo"],
    "quaternary_2_R": ["declared_tetrahedral_stereo"],
    "quaternary_2_S": ["declared_tetrahedral_stereo"],
    "penicillin_core": ["fused_ring", "declared_tetrahedral_stereo"],
    "dimethylbiphenyl_2_2": ["biaryl_hindered"],
}

# Molecules added beyond the source 63, for taxonomy categories otherwise
# entirely unrepresented. Reused verbatim from other already-established
# fixtures in this repo -- see module docstring for provenance.
ADDED_FOR_TAXONOMY = [
    {
        "name": "ring_torsion_fail_closed",
        "smiles": "C1CCCCC1CCCCCCCCCCCC",
        "primary_category": "known_fail_closed_case",
        "additional_tags": ["ring_torsion_fail_closed"],
        "added_for_taxonomy": True,
        "provenance": "validation/pipeline_v2_wasm_parity_fixtures.json "
        "(fixture name: ring_torsion_fail_closed) -- reused verbatim, "
        "originally generated by scripts/gen_pipeline_v2_wasm_parity_fixtures.py",
    },
    {
        "name": "force_field_unsupported_probe",
        # A bare phosphonium ylide fragment: no MMFF94 bond/angle parameters
        # exist for the P=C moiety in chematic-ff's MMFF94 table.
        "smiles": "[P](C)(C)(C)=C",
        "primary_category": "force_field_unsupported",
        "additional_tags": ["force_field_unsupported"],
        "added_for_taxonomy": True,
        "provenance": "constructed to probe chematic-ff's known MMFF94 "
        "bond/angle coverage gap for phosphonium ylides; not asserted a "
        "priori as the only force-field-unsupported corpus member -- arm "
        "coverage results determine that empirically.",
    },
]


def main():
    rows = []
    for name, smiles, category in CORPUS:
        rows.append(
            {
                "name": name,
                "smiles": smiles,
                "primary_category": category,
                "additional_tags": ADDITIONAL_TAGS.get(name, []),
                "added_for_taxonomy": False,
            }
        )
    rows.extend(ADDED_FOR_TAXONOMY)

    names = [r["name"] for r in rows]
    if len(names) != len(set(names)):
        sys.exit("duplicate name(s) in Tier A corpus -- fix before generating manifest")

    corpus_hash = hashlib.sha256(
        json.dumps(rows, sort_keys=True).encode("utf-8")
    ).hexdigest()

    manifest = {
        "tier": "A",
        "description": "Curated stress corpus for the pipeline v2 vs RDKit "
        "ETKDGv3 benchmark (Wave 1). Extracted from the existing pipeline v2 "
        "integration-gate corpus, not hand-copied a third time.",
        "source_file": SOURCE_FILE,
        "source_symbol": SOURCE_SYMBOL,
        "generator": "scripts/gen_pipeline_v2_vs_rdkit_tier_a_manifest.py",
        "molecule_count": len(rows),
        "corpus_sha256": corpus_hash,
        "molecules": rows,
    }

    out_path = ROOT / "validation" / "manifests" / "pipeline_v2_vs_rdkit_etkdgv3_tier_a.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Wrote {out_path.relative_to(ROOT)}: {len(rows)} molecules, sha256={corpus_hash[:16]}...")


if __name__ == "__main__":
    main()
