"""Compare NumPy destinations with RDKit list counters and fragment graphs."""

import hashlib
import itertools
import json
from pathlib import Path

import chematic as ch
import numpy as np
import rdkit
from rdkit import Chem


def main():
    cases = 0
    molecules = [
        "CC",
        "CCC",
        "CC(C)O",
        "CCOC",
        "CC(=O)O",
        "c1ccccc1CO",
        "[13CH3:7]CO",
        "CC[NH3+]",
        "CC.COC",
        "NCC(=O)O",
    ]
    for smiles in molecules:
        source = ch.from_smiles(smiles)
        graph = json.loads(source.to_moljson())
        params = Chem.SmilesParserParams()
        params.removeHs = False
        reference = Chem.MolFromSmiles(smiles, params)
        edges = [
            (int(b["source_id"][1:]) - 1, int(b["target_id"][1:]) - 1)
            for b in graph["bonds"]
        ]
        eligible = [
            i
            for i, (a, b) in enumerate(edges)
            if reference.GetBondBetweenAtoms(a, b).GetBondType() == Chem.BondType.SINGLE
            and not reference.GetBondBetweenAtoms(a, b).IsInRing()
        ]
        for ncuts in range(1, len(eligible) + 1):
            for selection in itertools.combinations(eligible, ncuts):
                for cuts in dict.fromkeys([selection, selection[::-1]]):
                    # Explicit equal dummy labels separate counter verification
                    # from native/SMILES bond endpoint orientation.
                    labels = [(17, 17)] * len(cuts)
                    rd_indices = [
                        reference.GetBondBetweenAtoms(*edges[i]).GetIdx() for i in cuts
                    ]
                    expected_counts = [99] * source.atom_count
                    expected = Chem.FragmentOnBonds(
                        reference,
                        rd_indices,
                        dummyLabels=labels,
                        cutsPerAtom=expected_counts,
                    )
                    for layout in ("plain", "stride", "reverse", "big_endian"):
                        storage = np.full(2 * source.atom_count, 77, dtype=np.uint32)
                        counts = {
                            "plain": storage[: source.atom_count],
                            "stride": storage[::2],
                            "reverse": storage[::2][::-1],
                            "big_endian": np.full(source.atom_count, 77, dtype=">u4"),
                        }[layout]
                        result = ch.fragment_on_bonds_checked(
                            source, cuts, dummy_labels=labels, cuts_per_atom=counts
                        )
                        assert (
                            counts.tolist()
                            == expected_counts
                            == list(result.cuts_per_atom)
                        )
                        if layout in ("stride", "reverse"):
                            assert np.all(storage[1::2] == 77)
                        parsed = Chem.MolFromSmiles(result.molecule.smiles, params)
                        assert parsed is not None
                        assert Chem.MolToSmiles(parsed) == Chem.MolToSmiles(expected)
                        cases += 1
    import chematic.chematic as native

    binary = Path(native.__file__)
    print(
        json.dumps(
            {
                "rdkit": rdkit.__version__,
                "numpy": np.__version__,
                "molecules": len(molecules),
                "nonempty_comparisons": cases,
                "mismatches": 0,
                "native_sha256": hashlib.sha256(binary.read_bytes()).hexdigest(),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
