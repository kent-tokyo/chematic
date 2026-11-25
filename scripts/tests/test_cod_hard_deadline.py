import multiprocessing
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'benchmarks'))
from cod_hard_deadline import run_engine


def test_hard_deadline_terminates_worker_and_next_call_succeeds():
    before = {child.pid for child in multiprocessing.active_children()}
    interrupted = run_engine('rdkit', 'C' * 300, 20, 0.2, 400)
    assert interrupted['status'] == 'timeout'
    assert interrupted['elapsed_with_cleanup'] < 3
    following = run_engine('rdkit', 'CCO', 1, 2.0, 400)
    assert following['status'] == 'ok'
    assert following['attempted'] == 1
    assert len(following['candidates']) == 1
    assert {child.pid for child in multiprocessing.active_children()} == before


def test_unsupported_chemical_does_not_prevent_rdkit():
    smiles = 'COc1ccc(C2=Cc3cc(Cl)nc[n+]3[B-](c3ccccc3)(c3ccccc3)O2)cc1'
    chematic = run_engine('chematic', smiles, 1, 2.0, 400)
    rdkit = run_engine('rdkit', smiles, 1, 2.0, 400)
    assert chematic['status'] == 'error'
    assert 'type assignment' in chematic['error']
    assert rdkit['status'] == 'ok'
    assert rdkit['candidates'][0][0] == -1
