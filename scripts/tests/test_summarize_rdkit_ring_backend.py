from scripts.summarize_rdkit_ring_backend import compare_triplet


def row(smarts):
    return {
        "input_index": 0,
        "smiles": "C1CC1",
        "status": "ok",
        "canonical": "C1CC1",
        "cip_atoms": {},
        "cip_bonds": {},
        "morgan_on_bits": [],
        "smarts": [smarts],
    }


def test_legacy_ring_probe_marks_reproduced_oracle_cell():
    old = row([[0]])
    new = row([])
    report = compare_triplet([old], [new], [row([[0]])], ["[R2]"])
    assert report["old_vs_new"]["smarts_changed_cells"] == 1
    assert report["legacy_reproduces_old_changed_cells"] == 1
    assert report["new_vs_new_legacy"]["smarts_changed_cells"] == 1
    assert report["old_vs_new_legacy"]["smarts_changed_cells"] == 0


def test_legacy_ring_probe_preserves_unexplained_change():
    report = compare_triplet([row([[0]])], [row([])], [row([[1]])], ["[R2]"])
    assert report["legacy_reproduces_old_changed_cells"] == 0
    assert report["old_vs_new_legacy"]["smarts_changed_cells"] == 1
