import copy
import json
import os
from pathlib import Path
import shutil
import sys

import pytest

from scripts.check_interchange_preservation import (
    ContractError, FIELD_NAMES, bounded_run, correspondence, evaluate,
    load_protocol, measure, publish_report, read_cjson, read_mol2, strict_json, summarize,
)

ROOT = Path(__file__).resolve().parents[2] / 'validation/interchange'
PROTOCOL = ROOT / 'protocol.json'


def sample():
    protocol = json.loads(PROTOCOL.read_bytes())
    case = protocol['cases'][0]
    candidate = read_mol2((ROOT / case['file']).read_text())
    return case['expected'], candidate


def check(expected, candidate, tolerance=1e-12, budget=10000):
    return evaluate(expected, candidate, coordinate_tolerance=tolerance,
                    charge_tolerance=tolerance, mapping_budget=budget)


def test_fixture_oracles_are_checked_before_measurement():
    protocol, cases, identity = load_protocol(PROTOCOL)
    assert len(cases) == 8
    assert len(identity['sha256']) == 64
    assert [c['id'] for c, _, _ in cases if not c['valid_in_declared_profile']] == [
        'nonfinite-charge', 'dangling-bond', 'unknown-bond']
    assert len(protocol['routes']) == 2


def test_mapping_handles_renumbering_and_reordering_without_fitting_tested_fields():
    expected, candidate = sample()
    # Candidate order: nitrogen, oxygen, methyl, carbonyl.
    permutation = [3, 2, 0, 1]
    inverse = {old: new for new, old in enumerate(permutation)}
    candidate['atoms'] = [candidate['atoms'][i] for i in permutation]
    for index, atom in enumerate(candidate['atoms']):
        atom['source_id'] = 100 + index
    candidate['bonds'] = [{'atoms': [inverse[i] for i in b['atoms']], 'order': b['order']}
                          for b in candidate['bonds']]
    result = check(expected, candidate)
    assert result['source_to_candidate'] == [2, 3, 1, 0]
    assert result['fields']['formal_charges'] == 'preserved'
    assert result['fields']['coordinates'] == 'preserved'
    assert result['fields']['source_atom_ids'] == 'lost'


@pytest.mark.parametrize('key,value,field', [
    ('formal_charge', 1, 'formal_charges'),
    ('partial_charge', -0.7, 'partial_charges'),
    ('coordinates', [5, 6, 7], 'coordinates'),
    ('name', 'changed', 'atom_names'),
    ('atom_type', 'C.2', 'atom_types'),
])
def test_mapped_field_changes_are_detected(key, value, field):
    expected, candidate = sample()
    candidate['atoms'][0][key] = value
    result = check(expected, candidate)
    assert result['correspondence'] == 'unique'
    assert result['fields'][field] == 'lost'


def test_bond_order_is_checked_after_unlabelled_topology_correspondence():
    expected, candidate = sample()
    candidate['bonds'][2]['order'] = 'aromatic'
    result = check(expected, candidate)
    assert result['fields']['connectivity'] == 'preserved'
    assert result['fields']['bond_orders'] == 'lost'


def test_atom_symmetry_is_not_resolved_using_coordinates_or_charges():
    expected, candidate = sample()
    expected['atoms'] = expected['atoms'][:2]
    expected['bonds'] = [{'atoms': [0, 1], 'order': 'single'}]
    candidate['atoms'] = candidate['atoms'][:2]
    candidate['bonds'] = expected['bonds']
    result = check(expected, candidate)
    assert result['correspondence'] == 'ambiguous'
    assert result['fields']['connectivity'] == 'preserved'
    assert result['fields']['coordinates'] == 'unverified'
    assert result['fields']['formal_charges'] == 'unverified'


def test_mapping_budget_exhaustion_is_explicit():
    expected, candidate = sample()
    mapping, reason = correspondence(expected, candidate, budget=1)
    assert mapping is None and reason == 'search_limit'
    assert check(expected, candidate, budget=1)['fields']['connectivity'] == 'unverified'


def test_precision_tolerance_is_separate_from_exact_f64_diagnostic():
    expected, candidate = sample()
    candidate['atoms'][0]['coordinates'][0] = 0.1235
    result = check(expected, candidate, tolerance=0.00005)
    assert result['fields']['coordinates'] == 'preserved'
    assert result['coordinate_bits_equal_diagnostic'] is False
    assert check(expected, candidate)['fields']['coordinates'] == 'lost'


def test_lost_connectivity_does_not_manufacture_correspondence():
    expected, candidate = sample()
    candidate['bonds'] = []
    result = check(expected, candidate)
    assert result['fields']['bond_count'] == 'lost'
    assert result['fields']['connectivity'] == 'lost'
    assert result['fields']['coordinates'] == 'unverified'
    assert result['source_to_candidate'] is None


def test_cjson_reader_distinguishes_partial_and_formal_charges():
    text = json.dumps({'atoms': {'elements': {'number': [7]}, 'formalCharges': [1],
                                'coords': {'3d': [0, 1, 2]}},
                       'partialCharges': {'USER_CHARGES': [0.8]},
                       'properties': {'name': 'ammonium projection'}})
    atom = read_cjson(text)['atoms'][0]
    assert atom['formal_charge'] == 1 and atom['partial_charge'] == 0.8
    assert read_cjson(text)['charge_method'] == 'USER_CHARGES'


@pytest.mark.parametrize('text', [
    '{"atoms":{"elements":{"number":[6]},"formalCharges":[]}}',
    '{"atoms":{"elements":{"number":[6]},"coords":{"3d":[0,0]}}}',
    '{"atoms":{"elements":{"number":[6]}},"bonds":{"connections":{"index":[0,1]},"order":[1]}}',
    '{"atoms":{"elements":{"number":[6]}},"bonds":{"connections":{"index":[0,0]},"order":[1]}}',
    '{"atoms":{"elements":{"number":[6]}},"atoms":{}}',
    '{"atoms":{"elements":{"number":[true]}}}',
])
def test_independent_reader_rejects_malformed_output(text):
    with pytest.raises((ValueError, KeyError, TypeError)):
        read_cjson(text)


@pytest.mark.parametrize('text', ['{"x":NaN}', '{"x":1,"x":2}'])
def test_json_ambiguity_is_rejected(text):
    with pytest.raises(ContractError):
        strict_json(text)


def test_misdeclared_oracle_is_rejected(tmp_path):
    shutil.copytree(ROOT, tmp_path / 'corpus')
    path = tmp_path / 'corpus/protocol.json'
    protocol = json.loads(path.read_bytes())
    protocol['cases'][0]['expected']['atoms'][0]['partial_charge'] = -0.7
    path.write_text(json.dumps(protocol))
    with pytest.raises(ContractError, match='fixture disagrees'):
        load_protocol(path)


@pytest.mark.parametrize('change', ['duplicate_id', 'escape', 'bad_tolerance', 'budget'])
def test_invalid_protocol_is_rejected(tmp_path, change):
    shutil.copytree(ROOT, tmp_path / 'corpus')
    path = tmp_path / 'corpus/protocol.json'
    protocol = json.loads(path.read_bytes())
    if change == 'duplicate_id':
        protocol['cases'][1]['id'] = protocol['cases'][0]['id']
    elif change == 'escape':
        protocol['cases'][0]['file'] = '../../outside.mol2'
    elif change == 'bad_tolerance':
        protocol['coordinate_abs_tolerance']['mol2'] = -1
    else:
        protocol['limits']['max_atoms'] = 10000
    path.write_text(json.dumps(protocol))
    with pytest.raises(ContractError):
        load_protocol(path)


def test_output_limit_is_enforced_during_capture():
    result = bounded_run([sys.executable, '-c', "import os; os.write(1,b'x'*100000); os.write(2,b'y'*100000)"], max_output=1000)
    assert result['resource_status'] == 'output_limit'
    assert len(result['stdout']) + len(result['stderr']) == 1000


def test_timeout_includes_blocked_stdin_and_descendant_pipes():
    result = bounded_run([sys.executable, '-c', 'import time; time.sleep(10)'], b'x' * 1048576, timeout=0.2)
    assert result['resource_status'] == 'timeout'
    code = "import subprocess,sys; subprocess.Popen([sys.executable,'-c','import time; time.sleep(10)'])"
    result = bounded_run([sys.executable, '-c', code], timeout=0.2)
    assert result['resource_status'] == 'timeout'


def test_nonzero_exit_and_both_pipes_are_retained():
    result = bounded_run([sys.executable, '-c', "import sys; print('out'); print('err',file=sys.stderr); sys.exit(7)"])
    assert result['resource_status'] == 'completed'
    assert result['exit_code'] == 7
    assert result['stdout'] == b'out\n' and result['stderr'] == b'err\n'


def test_safe_refusal_does_not_count_as_preserved_coverage(tmp_path):
    protocol, cases, _ = load_protocol(PROTOCOL)
    case, source, _ = cases[0]
    report = {'schema': 'chematic.conversion-report.v1', 'status': 'refused', 'output': None}
    command = [sys.executable, '-c', f'import sys; print({json.dumps(report)!r}); sys.exit(2)']
    row = measure('chematic', command, source, case, protocol['routes'][0], protocol, tmp_path)
    assert row['outcome'] == 'safe_refusal'
    summary = summarize([row])[0]
    assert summary['input_count'] == summary['terminal_count'] == 1
    assert summary['outcomes'] == {'safe_refusal': 1}
    assert set(summary['fields']) == set(FIELD_NAMES)
    assert all(v == {'no_candidate': 1} for v in summary['fields'].values())


def test_atomic_publication_never_overwrites_existing_record(tmp_path):
    path = tmp_path / 'report.json'
    publish_report(path, {'first': True})
    with pytest.raises(FileExistsError):
        publish_report(path, {'second': True})
    assert json.loads(path.read_bytes()) == {'first': True}
    assert list(tmp_path.iterdir()) == [path]
    with pytest.raises(ValueError):
        publish_report(tmp_path / 'invalid.json', {'bad': float('nan')})
    assert list(tmp_path.iterdir()) == [path]

@pytest.mark.parametrize('mismatch', ['schematic', 'openbabel'])
def test_engine_artifact_mismatch_is_rejected_before_creating_outputs(tmp_path, mismatch):
    import hashlib
    import subprocess
    binary = Path(sys.executable).resolve()
    data = binary.read_bytes()
    digest = hashlib.sha256(data).hexdigest()
    evidence = tmp_path / 'build.json'
    evidence.write_text(json.dumps({'schema': 'chematic.development-verification.v1',
                                    'workspace_version': '1.0.20',
                                    'artifacts': [{'sha256': 'bad' if mismatch == 'schematic' else digest,
                                                   'size_bytes': len(data)}]}))
    output = tmp_path / 'report.json'
    command = [sys.executable, str(ROOT.parents[1] / 'scripts/check_interchange_preservation.py'),
               '--protocol', str(PROTOCOL), '--chematic', str(binary),
               '--chematic-evidence', str(evidence), '--obabel', str(binary),
               '--expected-obabel-version', '3.2.1',
               '--expected-obabel-sha256', 'bad' if mismatch == 'openbabel' else digest,
               '--output', str(output)]
    result = subprocess.run(command, capture_output=True, timeout=10)
    assert result.returncode == 2
    assert b'differs from pinned SHA-256' in result.stderr if mismatch == 'openbabel' else b'does not match selected build evidence' in result.stderr
    assert not output.exists()
    assert not output.with_suffix('.json.artifacts').exists()
