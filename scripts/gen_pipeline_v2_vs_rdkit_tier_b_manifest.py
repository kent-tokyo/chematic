#!/usr/bin/env python3
"""Generate the Tier B (fixed drug-like) corpus manifest for the pipeline v2
vs RDKit ETKDGv3 benchmark (Wave 1 of the "RDKit alternative" program).

Source: `scripts/descriptor_census_corpus.smi` -- a 5,000-molecule ChEMBL
corpus already committed to this repo for reproducibility (downloaded via
`scripts/download_chembl_smiles.py --count 5000`, documented in
`docs/descriptor_census_rfc.md`). Deliberately reused rather than
re-downloaded: it is already committed (no `~/Downloads` dependency), already
license-cleared, and its own provenance is already on record. This script
draws its OWN, independently-defined 200-molecule subset from it (different
selection rule and different resulting set than the descriptor census's own
use of the full 5,000) -- not a copy of any other benchmark's Tier B.

License: ChEMBL database, CC Attribution-ShareAlike 3.0
(https://www.ebi.ac.uk/about/terms-of-use) -- canonical SMILES only
redistributed, no other ChEMBL fields.

Selection rule (fully deterministic, no RNG):
  1. Read `scripts/descriptor_census_corpus.smi` in file order.
  2. Parse each line with RDKit; drop parse failures (counted, not silently
     skipped -- see `parse_failure_count` in the output manifest).
  3. Drop multi-fragment entries (a "." in the SMILES) -- salts/disconnected
     structures are out of scope for a single-conformer 3D embedding
     benchmark.
  4. Element filter: heavy atoms restricted to
     {C, N, O, F, P, S, Cl, Br, I} -- standard organic drug-like scope,
     documented explicitly rather than left implicit.
  5. Heavy-atom range: [8, 60] inclusive -- excludes both trivial fragments
     and outliers far outside typical drug-like size.
  6. Dedup by RDKit canonical SMILES (first occurrence wins; later duplicates
     dropped and counted).
  7. Take the first `TARGET_COUNT` molecules surviving steps 2-6, in file
     order -- deterministic, no random sampling.

Run: `.venv/bin/python scripts/gen_pipeline_v2_vs_rdkit_tier_b_manifest.py`
Output: `validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_b.json`
"""

import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).parent.parent
SOURCE_PATH = ROOT / "scripts" / "descriptor_census_corpus.smi"
OUT_PATH = ROOT / "validation" / "manifests" / "pipeline_v2_vs_rdkit_etkdgv3_tier_b.json"

TARGET_COUNT = 200
ALLOWED_ELEMENTS = {"C", "N", "O", "F", "P", "S", "Cl", "Br", "I"}
MIN_HEAVY_ATOMS = 8
MAX_HEAVY_ATOMS = 60


def main():
    try:
        from rdkit import Chem
    except ImportError:
        sys.exit("rdkit not installed. pip install rdkit")

    source_bytes = SOURCE_PATH.read_bytes()
    source_sha256 = hashlib.sha256(source_bytes).hexdigest()

    lines = source_bytes.decode("utf-8").splitlines()

    selected = []
    seen_canonical = set()
    parse_failures = 0
    multi_fragment_dropped = 0
    element_filtered = 0
    size_filtered = 0
    dedup_dropped = 0
    scanned = 0

    for line in lines:
        smiles = line.strip()
        if not smiles:
            continue
        scanned += 1
        if len(selected) >= TARGET_COUNT:
            break  # deterministic: stop as soon as target is met, don't over-scan

        if "." in smiles:
            multi_fragment_dropped += 1
            continue

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            parse_failures += 1
            continue

        elements = {atom.GetSymbol() for atom in mol.GetAtoms()}
        if not elements.issubset(ALLOWED_ELEMENTS):
            element_filtered += 1
            continue

        heavy_count = mol.GetNumAtoms()
        if not (MIN_HEAVY_ATOMS <= heavy_count <= MAX_HEAVY_ATOMS):
            size_filtered += 1
            continue

        canonical = Chem.MolToSmiles(mol)
        if canonical in seen_canonical:
            dedup_dropped += 1
            continue
        seen_canonical.add(canonical)

        selected.append(
            {
                "name": f"chembl_tier_b_{len(selected):04d}",
                "smiles": canonical,
                "heavy_atom_count": heavy_count,
                "primary_category": "drug_like",
            }
        )

    if len(selected) < TARGET_COUNT:
        print(
            f"WARNING: only {len(selected)}/{TARGET_COUNT} molecules survived "
            f"filtering after scanning all {scanned} source lines -- using "
            f"{len(selected)} as the actual Tier B corpus size.",
            file=sys.stderr,
        )

    corpus_hash = hashlib.sha256(
        json.dumps(selected, sort_keys=True).encode("utf-8")
    ).hexdigest()

    manifest = {
        "tier": "B",
        "description": "Fixed drug-like corpus for the pipeline v2 vs RDKit "
        "ETKDGv3 benchmark (Wave 1), deterministically subsetted from an "
        "already-committed ChEMBL corpus.",
        "source_file": "scripts/descriptor_census_corpus.smi",
        "source_sha256": source_sha256,
        "source_license": "ChEMBL database, CC BY-SA 3.0 "
        "(https://www.ebi.ac.uk/about/terms-of-use); canonical SMILES only.",
        "source_acquisition": "scripts/download_chembl_smiles.py --count 5000 "
        "(already committed to this repo; not re-downloaded for this "
        "benchmark -- see docs/descriptor_census_rfc.md)",
        "generator": "scripts/gen_pipeline_v2_vs_rdkit_tier_b_manifest.py",
        "target_count": TARGET_COUNT,
        "selection_rule": [
            "read source file in line order",
            "drop parse failures (RDKit)",
            "drop multi-fragment SMILES ('.')",
            f"element filter: heavy atoms in {sorted(ALLOWED_ELEMENTS)}",
            f"heavy-atom range: [{MIN_HEAVY_ATOMS}, {MAX_HEAVY_ATOMS}]",
            "dedup by RDKit canonical SMILES (first occurrence wins)",
            f"take first {TARGET_COUNT} survivors in file order (deterministic, no RNG)",
        ],
        "source_lines_scanned": scanned,
        "parse_failure_count": parse_failures,
        "multi_fragment_dropped_count": multi_fragment_dropped,
        "element_filtered_count": element_filtered,
        "size_filtered_count": size_filtered,
        "dedup_dropped_count": dedup_dropped,
        "molecule_count": len(selected),
        "corpus_sha256": corpus_hash,
        "molecules": selected,
    }

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUT_PATH.write_text(json.dumps(manifest, indent=2) + "\n")
    print(
        f"Wrote {OUT_PATH.relative_to(ROOT)}: {len(selected)} molecules "
        f"(scanned {scanned} source lines), corpus_sha256={corpus_hash[:16]}..."
    )


if __name__ == "__main__":
    main()
