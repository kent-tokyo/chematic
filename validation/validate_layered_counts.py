#!/usr/bin/env python3
"""Live RDKit reference comparisons for typed caller-owned Layered counters."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import rdkit
from rdkit import Chem, DataStructs
from chematic import rdkit_compat as Ch


SMILES = [
    "CCO",
    "CC(C)O",
    "c1ccccc1",
    "C1CCCCC1",
    "c1ccncc1",
    "CC(=O)N",
    "C1CC2CCC1C2",
    "C1COC1",
    "C1CCC2(CC1)CCCC2",
    "[Na+].CC(=O)[O-]",
    "C#N",
    "[13CH3]O",
    "C1=CC=CC=C1Cl",
    "CCN(C)C",
    "C[C@H](O)F",
    "C",
]


def main():
    comparisons = 0
    view_checks = 0
    accumulated_checks = 0
    digest = hashlib.sha256()
    for smiles in SMILES:
        mol = Ch.MolFromSmiles(smiles)
        reference = Chem.MolFromSmiles(smiles)
        assert mol.GetNumAtoms() == reference.GetNumAtoms()
        n = mol.GetNumAtoms()
        for width in [1, 9, 257, 2048]:
            for flag in [0, 1, 4, 0x3F]:
                for roots, branched in [
                    (None, True),
                    ([0], True),
                    ([0, 0], True),
                    ([0], False),
                ]:
                    for masked in [False, True]:
                        mask = Ch.ExplicitBitVect(width) if masked else None
                        ref_mask = (
                            DataStructs.ExplicitBitVect(width) if masked else None
                        )
                        if masked:
                            for bit in range(0, width, 3):
                                mask.SetBit(bit)
                                ref_mask.SetBit(bit)
                        seed = [2**32 - 2] + [7] * (n - 1) + [91, 92]
                        counts = seed.copy()
                        ref_counts = seed.copy()
                        kwargs = {
                            "layerFlags": flag,
                            "minPath": 1,
                            "maxPath": 4,
                            "fpSize": width,
                            "branchedPaths": branched,
                        }
                        if roots is not None:
                            kwargs["fromAtoms"] = roots
                        expected = Chem.LayeredFingerprint(
                            reference,
                            atomCounts=ref_counts,
                            setOnlyBits=ref_mask,
                            **kwargs,
                        )
                        actual = Ch.LayeredFingerprint(
                            mol, atomCounts=counts, setOnlyBits=mask, **kwargs
                        )
                        assert actual.ToBitString() == expected.ToBitString(), (
                            smiles,
                            kwargs,
                            masked,
                            "bits",
                        )
                        assert counts == ref_counts, (
                            smiles,
                            kwargs,
                            masked,
                            "counts",
                            counts,
                            ref_counts,
                        )
                        # A second invocation must accumulate on the same public destination.
                        again = Chem.LayeredFingerprint(
                            reference,
                            atomCounts=ref_counts,
                            setOnlyBits=ref_mask,
                            **kwargs,
                        )
                        Ch.LayeredFingerprint(
                            mol, atomCounts=counts, setOnlyBits=mask, **kwargs
                        )
                        assert counts == ref_counts
                        accumulated_checks += 1
                        # Native-endian strided and non-native-endian reversed uint32 views.
                        for endian, reverse in [("<u4", False), (">u4", True)]:
                            storage = np.full(2 * (n + 2), 123, dtype=endian)
                            view = storage[::-2] if reverse else storage[::2]
                            view[:] = seed
                            Ch.LayeredFingerprint(
                                mol, atomCounts=view, setOnlyBits=mask, **kwargs
                            )
                            assert view.tolist() == [
                                int(v)
                                for v in Chem_counts_once(
                                    reference, ref_mask, seed, kwargs
                                )
                            ]
                            assert (
                                storage[::2].tolist() == [123] * (n + 2)
                                if reverse
                                else storage[1::2].tolist() == [123] * (n + 2)
                            )
                            view_checks += 1
                        digest.update(
                            json.dumps(
                                [
                                    smiles,
                                    width,
                                    flag,
                                    roots,
                                    branched,
                                    masked,
                                    counts,
                                    actual.ToBitString(),
                                ]
                            ).encode()
                        )
                        comparisons += 1
    print(
        json.dumps(
            {
                "rdkit_version": rdkit.__version__,
                "numpy_version": np.__version__,
                "molecules": len(SMILES),
                "reference_comparisons": comparisons,
                "second_call_accumulation_checks": accumulated_checks,
                "numpy_view_checks": view_checks,
                "corpus_sha256": digest.hexdigest(),
            },
            indent=2,
            sort_keys=True,
        )
    )


def Chem_counts_once(mol, mask, seed, kwargs):
    counts = seed.copy()
    Chem.LayeredFingerprint(mol, atomCounts=counts, setOnlyBits=mask, **kwargs)
    return counts


if __name__ == "__main__":
    main()
