"""Cross-language shared-corpus test for the promoted RDKit-exact Morgan/ECFP API
(`Mol.rdkit_ecfp4` / `Mol.rdkit_ecfp4_detail` / `Mol.rdkit_ecfp_config` /
`Mol.rdkit_ecfp_config_detail`).

Reads the SAME `validation/ecfp4_rdkit_stable_api_fixtures.json` file the Rust test
(`crates/chematic-fp/tests/rdkit_morgan_stable_api_fixtures.rs`) and the WASM test
(`crates/chematic-wasm/tests/rdkit_ecfp4_stable_api.test.mjs`) read -- generated once
from a live RDKit oracle by `scripts/gen_ecfp4_rdkit_stable_api_fixtures.py` (RDKit
version/commit recorded in the file itself). Does not import RDKit itself -- the
oracle values are already baked into the JSON, so this test has no RDKit dependency
and always runs (unlike the `pytest.importorskip("rdkit")`-gated differential tests
elsewhere in this suite).
"""
import json
import os

import pytest

import chematic

_CORPUS_PATH = os.path.join(
    os.path.dirname(__file__), "..", "..", "..",
    "validation", "ecfp4_rdkit_stable_api_fixtures.json",
)

with open(_CORPUS_PATH) as f:
    _CORPUS = json.load(f)

_FIXTURES = _CORPUS["fixtures"]


def _fingerprint_bytes_to_bit_list(fp_bytes):
    bits = []
    for byte_idx, byte in enumerate(fp_bytes):
        for bit in range(8):
            if byte & (1 << bit):
                bits.append(byte_idx * 8 + bit)
    return sorted(bits)


def test_corpus_has_real_coverage():
    assert _CORPUS["rdkit_version"], "corpus must record the RDKit version"
    assert len(_FIXTURES) >= 30
    ok = [fx for fx in _FIXTURES if fx["expect"] == "ok"]
    err = [fx for fx in _FIXTURES if fx["expect"] == "error"]
    assert len(ok) >= 30
    assert len(err) >= 1


@pytest.mark.parametrize("fx", [fx for fx in _FIXTURES if fx["expect"] == "ok"], ids=lambda fx: fx["id"])
def test_rdkit_ecfp4_matches_shared_oracle_corpus(fx):
    mol = chematic.from_smiles(fx["smiles"])

    fp_bytes = mol.rdkit_ecfp4()
    assert isinstance(fp_bytes, bytes)
    assert len(fp_bytes) == 256
    assert _fingerprint_bytes_to_bit_list(fp_bytes) == fx["folded_bits"]

    fp2, sparse_counts, raw_bit_info, folded_bit_info = mol.rdkit_ecfp4_detail()
    assert fp2 == fp_bytes

    expected_sparse = {int(k): v for k, v in fx["sparse_counts"].items()}
    assert sparse_counts == expected_sparse

    expected_raw = {int(k): sorted(tuple(p) for p in v) for k, v in fx["raw_bit_info"].items()}
    got_raw = {k: sorted(tuple(p) for p in v) for k, v in raw_bit_info.items()}
    assert got_raw == expected_raw

    expected_folded = {int(k): sorted(tuple(p) for p in v) for k, v in fx["folded_bit_info"].items()}
    got_folded = {k: sorted(tuple(p) for p in v) for k, v in folded_bit_info.items()}
    assert got_folded == expected_folded


@pytest.mark.parametrize("fx", [fx for fx in _FIXTURES if fx["expect"] == "error"], ids=lambda fx: fx["id"])
def test_rdkit_ecfp4_raises_value_error_on_preprocessing_failure(fx):
    mol = chematic.from_smiles(fx["smiles"])
    with pytest.raises(ValueError):
        mol.rdkit_ecfp4()
    with pytest.raises(ValueError):
        mol.rdkit_ecfp4_detail()
    with pytest.raises(ValueError):
        mol.rdkit_ecfp_config(radius=2, nbits=2048)


@pytest.mark.parametrize("fx", [fx for fx in _FIXTURES if fx["expect"] == "ok"], ids=lambda fx: fx["id"])
def test_rdkit_ecfp_config_radius_axis_matches_shared_oracle_corpus(fx):
    mol = chematic.from_smiles(fx["smiles"])
    for cell in fx["radius_axis"]:
        fp_bytes = mol.rdkit_ecfp_config(radius=cell["radius"], nbits=2048)
        assert _fingerprint_bytes_to_bit_list(fp_bytes) == cell["folded_bits"], (
            f"fixture {fx['id']} radius={cell['radius']}"
        )


@pytest.mark.parametrize("fx", [fx for fx in _FIXTURES if fx["expect"] == "ok"], ids=lambda fx: fx["id"])
def test_rdkit_ecfp_config_fp_size_axis_matches_shared_oracle_corpus(fx):
    mol = chematic.from_smiles(fx["smiles"])
    for cell in fx["fp_size_axis"]:
        fp_bytes = mol.rdkit_ecfp_config(radius=2, nbits=cell["fp_size"])
        assert len(fp_bytes) == cell["fp_size"] // 8
        assert _fingerprint_bytes_to_bit_list(fp_bytes) == cell["folded_bits"], (
            f"fixture {fx['id']} fp_size={cell['fp_size']}"
        )


def test_rdkit_ecfp_config_default_matches_rdkit_ecfp4():
    mol = chematic.from_smiles("c1ccc2ccccc2c1")  # naphthalene
    assert mol.rdkit_ecfp_config() == mol.rdkit_ecfp4()


def test_rdkit_ecfp_config_rejects_unsupported_radius_and_nbits():
    mol = chematic.from_smiles("c1ccccc1")
    with pytest.raises(ValueError):
        mol.rdkit_ecfp_config(radius=4, nbits=2048)
    with pytest.raises(ValueError):
        mol.rdkit_ecfp_config(radius=2, nbits=64)


def test_rdkit_ecfp4_differs_from_legacy_ecfp4():
    # Different hash function/byte layout by design -- never silently interchanged.
    mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")  # aspirin
    assert mol.rdkit_ecfp4() != mol.ecfp4()
