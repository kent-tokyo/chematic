"""Pin default Pattern bits and actual no-op option behavior against RDKit."""

import hashlib
import json
from pathlib import Path

import chematic.chematic as native
from chematic import rdkit_compat as C
import numpy as np
import rdkit
from rdkit import Chem, DataStructs


def main():
    paths = [
        Path("scripts/descriptor_census_corpus.smi"),
        Path("scripts/chembl_accuracy_corpus_4999.smi"),
    ]
    samples = [
        line.split()[0] for path in paths for line in path.read_text().splitlines()[:32]
    ]
    samples += [
        "C",
        "CCO",
        "CC(=O)O",
        "c1ccccc1",
        "c1ccncc1",
        "CC[NH3+]",
        "[13CH3:7]CO",
        "C/C=C/C",
        "N[C@@H](C)O",
        "CC.O",
        "C#N",
        "[Na+].[Cl-]",
    ]
    comparisons = 0
    option_checks = 0
    mismatches = []
    for source_index, smiles in enumerate(samples):
        reference = Chem.MolFromSmiles(smiles)
        candidate = C.MolFromSmiles(smiles)
        assert reference is not None and candidate is not None
        for width in (1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048):
            rd_seed = [77] * (reference.GetNumAtoms() + 2)
            rd_mask = DataStructs.ExplicitBitVect(width)
            # An empty mask still produces normal bits upstream: no filtering.
            expected = Chem.PatternFingerprint(reference, width, rd_seed, rd_mask)
            assert rd_seed == [77] * len(rd_seed)
            expected_bits = tuple(expected.GetOnBits())
            seeds = [77] * (candidate.GetNumAtoms() + 2)
            storage = np.full(2 * len(seeds), 77, dtype=np.uint32)
            mask = C.ExplicitBitVect(width)
            for counts, supplied_mask in (
                (None, None),
                (seeds, mask),
                (storage[::-2], mask),
            ):
                observed = C.PatternFingerprint(
                    candidate, np.int64(width), counts, supplied_mask
                )
                actual_bits = tuple(observed.GetOnBits())
                if actual_bits != expected_bits:
                    mismatches.append(
                        {
                            "sample": source_index,
                            "smiles": smiles,
                            "width": width,
                            "counter": type(counts).__name__,
                        }
                    )
                assert seeds == [77] * len(seeds)
                assert np.all(storage == 77)
                comparisons += 1
                option_checks += int(counts is not None)
    summary = {
        "rdkit": rdkit.__version__,
        "numpy": np.__version__,
        "samples": len(samples),
        "widths": 12,
        "input_profiles": 3,
        "bit_vector_comparisons": comparisons,
        "no_op_option_checks": option_checks,
        "mismatch_count": len(mismatches),
        "mismatches": mismatches,
        "corpus_sha256": {
            str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths
        },
        "native_sha256": hashlib.sha256(Path(native.__file__).read_bytes()).hexdigest(),
    }
    print(json.dumps(summary, indent=2))
    if mismatches:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
