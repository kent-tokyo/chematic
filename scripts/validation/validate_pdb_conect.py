#!/usr/bin/env python3
"""Compare explicit indexed PDB graphs to pinned RDKit, without perception.

Requires the newly built Python module and Node WASM module. The optional
--issue-pdb file is the unedited PDB block from RDKit #9581. Raw inputs and
all compared outputs are retained in a deterministic gzip packet.
"""
import argparse
import gzip
import hashlib
import itertools
import json
import platform
import random
import subprocess
from pathlib import Path

import chematic
import chematic.chematic as extension
from rdkit import Chem, RDLogger, rdBase

ROOT = Path(__file__).resolve().parents[2]
FIXTURE = ROOT / 'validation/rdkit_issues/pdb-conect/sulfur.pdb'


def sha(data):
    return hashlib.sha256(data).hexdigest()


def graph(mol):
    assert mol is not None
    return {
        'atoms': [{'atomic_number': a.GetAtomicNum(), 'charge': a.GetFormalCharge()} for a in mol.GetAtoms()],
        'bonds': sorted([min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()), max(b.GetBeginAtomIdx(), b.GetEndAtomIdx()), int(b.GetBondTypeAsDouble())] for b in mol.GetBonds()),
    }


def reference(text):
    mol = Chem.MolFromPDBBlock(text, sanitize=False, removeHs=False, proximityBonding=False)
    result = graph(mol)
    conf = mol.GetConformer()
    result['coords'] = [[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z] for i in range(mol.GetNumAtoms())]
    result['serials'] = [a.GetPDBResidueInfo().GetSerialNumber() for a in mol.GetAtoms()]
    return result


def native(text, extended):
    f = chematic.from_pdb_conect(text, allow_extended_targets=extended)
    # A graph-only writer supplies indexed atoms/charges/edges without invoking
    # compatibility preparation or removing explicit H. Coordinates are checked
    # separately from the exact source frame.
    block = f['mol'].to_mol_v3000([[0., 0.]] * f['mol'].atom_count)
    g = graph(Chem.MolFromMolBlock(block, sanitize=False, removeHs=False))
    return dict(g, coords=f['coords'], serials=f['serials'])


def atom(serial, symbol, coord, charge=0):
    suffix = f'{abs(charge)}{"+" if charge > 0 else "-"}' if charge else '  '
    return f'HETATM{serial:5d}  CA  TST A   1    ' + ''.join(f'{x:8.3f}' for x in coord) + f'  1.00  0.00          {symbol:>2s}{suffix}'


def corpus():
    base = FIXTURE.read_text(encoding='utf-8')
    rows = [s for s in base.splitlines() if s.startswith('HETATM')]
    edges = [(10, 20, 2), (10, 30, 2), (10, 40, 1), (10, 50, 1)]
    for seed, spread, reciprocal, reverse, crlf in itertools.product(range(4), (False, True), (False, True), (False, True), (False, True)):
        order = list(rows)
        random.Random(seed).shuffle(order)
        if spread:
            order = [s[:30] + ''.join(f'{float(s[i:i+8])*10+9:8.3f}' for i in (30,38,46)) + s[54:] for s in order]
        links = []
        for a, b, n in edges:
            links.append('CONECT' + f'{a:5d}' + f'{b:5d}' * n)
            if reciprocal:
                links.append('CONECT' + f'{b:5d}' + f'{a:5d}' * n)
        if reverse:
            links.reverse()
        text = '\n'.join(order + links + ['END']) + '\n'
        if crlf:
            text = text.replace('\n', '\r\n')
        yield {'id': f'sulfur-{seed}-{spread}-{reciprocal}-{reverse}-{crlf}', 'input': text, 'extended': False}
    for n, reciprocal, split, symbols in itertools.product(range(1, 4), (False, True), (False, True), (('C','N'), ('Ca','Cl'), ('Br','H'))):
        # This is an indexed-graph contract, not a claim these artificial ion
        # and multiplicity combinations have chemically valid valence.
        rows = [atom(37, symbols[0], (0, 0, 0), 1), atom(12, symbols[1], (8, -2, 1), -1)]
        links = ['CONECT   37' + '   12' * n] if not split else ['CONECT   37   12'] * n
        if reciprocal:
            links.append('CONECT   12' + '   37' * n)
        yield {'id': f'pair-{n}-{reciprocal}-{split}-{symbols}', 'input': '\n'.join(rows + links + ['END']) + '\n', 'extended': False}
    yield {'id': 'no-edges', 'input': '\n'.join(rows) + '\nEND\n', 'extended': False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--wasm-module', type=Path, required=True)
    parser.add_argument('--issue-pdb', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert rdBase.rdkitVersion == '2026.09.1', rdBase.rdkitVersion
    RDLogger.DisableLog('rdApp.*')
    cases = list(corpus())
    issue = None
    if args.issue_pdb:
        data = args.issue_pdb.read_bytes()
        text = data.decode('utf-8')
        assert Chem.MolFromPDBBlock(text, sanitize=True, removeHs=False, proximityBonding=True) is None
        checked = Chem.MolFromPDBBlock(text, sanitize=True, removeHs=False, proximityBonding=False)
        assert checked is not None
        sulfur = checked.GetAtomWithIdx(0)
        assert sulfur.GetAtomicNum() == 16 and sum(b.GetBondTypeAsDouble() for b in sulfur.GetBonds()) == 6
        try:
            chematic.from_pdb_conect(text)
        except ValueError as error:
            assert 'extended CONECT' in str(error)
        else:
            raise AssertionError('legacy/extended targets must require explicit opt-in')
        legacy, _ = chematic.from_pdb(text)
        legacy_graph = graph(Chem.MolFromMolBlock(legacy.to_mol_v3000([[0,0]] * legacy.atom_count), sanitize=False, removeHs=False))
        legacy_sulfur = [edge for edge in legacy_graph['bonds'] if 0 in edge[:2]]
        cases.append({'id':'rdkit-9581-original', 'input':text, 'extended':True, 'legacy_native_graph':legacy_graph})
        issue = {'url':'https://github.com/rdkit/rdkit/issues/9581', 'input_sha256':sha(data), 'default_proximity_sanitized_result':None, 'proximity_disabled_sulfur_bond_order_sum':6, 'native_default':'extended targets refused', 'native_opt_in':'exact indexed graph comparison', 'legacy_native': {'atoms':len(legacy_graph['atoms']), 'bonds':len(legacy_graph['bonds']), 'sulfur_edges':legacy_sulfur, 'sulfur_bond_order_sum':sum(e[2] for e in legacy_sulfur)}}
    node = subprocess.run(['node', str(ROOT/'scripts/validation/validate_pdb_conect_node.mjs'), str(args.wasm_module.resolve())], input=json.dumps(cases), text=True, capture_output=True, check=True)
    node_rows = json.loads(node.stdout)
    assert len(node_rows) == len(cases)
    packet = []
    for case, js in zip(cases, node_rows):
        want = reference(case['input'])
        got = native(case['input'], case['extended'])
        assert got == want, (case['id'], got, want)
        js['frame']['bonds'] = sorted([min(a,b), max(a,b), n] for a,b,n in js['frame']['bonds'])
        assert js['frame'] == want, (case['id'], js, want)
        handle = graph(Chem.MolFromMolBlock(js['mol_block'], sanitize=False, removeHs=False))
        assert handle == {k:want[k] for k in ('atoms','bonds')}, case['id']
        packet.append(dict(case, rdkit=want, python=got, wasm=js))
    args.output.mkdir(parents=True, exist_ok=True)
    raw = gzip.compress(json.dumps(packet, sort_keys=True, separators=(',',':')).encode(), mtime=0)
    (args.output/'comparisons.json.gz').write_bytes(raw)
    sources = ['crates/chematic-mol/src/pdb_conect.rs', 'crates/chematic-py/src/formats/pdb_conect.rs', 'crates/chematic-wasm/src/pdb_conect.rs', 'scripts/validation/validate_pdb_conect.py', 'scripts/validation/validate_pdb_conect_node.mjs', 'validation/rdkit_issues/pdb-conect/sulfur.pdb']
    report = {
        'schema': 'chematic.pdb-conect-validation.v1',
        'scope': 'explicit indexed atoms/charges/bond orders/coordinates/serials; no sanitizer, aromaticity, implicit-H, stereo, general PDB or speed claim',
        'source_base': subprocess.check_output(['git','rev-parse','HEAD'], cwd=ROOT, text=True).strip(),
        'source_sha256': {s:sha((ROOT/s).read_bytes()) for s in sources},
        'runtime': {'python':platform.python_version(), 'rdkit':rdBase.rdkitVersion, 'platform':platform.platform(), 'python_extension_sha256':sha(Path(extension.__file__).read_bytes()), 'wasm_sha256':sha(args.wasm_module.with_name('chematic_wasm_bg.wasm').read_bytes())},
        'comparisons': {'cases':len(cases), 'python_passed':len(cases), 'wasm_frame_passed':len(cases), 'wasm_handle_graph_passed':len(cases), 'failures':0},
        'issue':issue,
        'packet': {'path':'comparisons.json.gz', 'sha256':sha(raw)},
    }
    (args.output/'summary.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
