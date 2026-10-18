"""Cross-read freshly compiled Rust output; run with RDKit on PYTHONPATH."""
from pathlib import Path
import io
import json
import subprocess
import rdkit
from rdkit import Chem, RDLogger

RDLogger.DisableLog('rdApp.warning')
root = Path(__file__).resolve().parent
records = [json.loads(line.split('TITLE_ORACLE ', 1)[1])
           for line in (root / 'after.log').read_text().splitlines()
           if line.startswith('TITLE_ORACLE ')]
assert len(records) == 126, len(records)
base_bodies = {r['writer']: r['block'].split('\n', 1)[1]
               for r in records if r['source_title'] == 'plain'}
summary = {'rdkit_version': rdkit.__version__, 'cases': 0, 'sdf_property_reads': 0,
           'unchanged_non_title_bodies': 0, 'mismatches': []}
for record in records:
    block = record['block']
    first, body = block.split('\n', 1)
    assert first == record['expected_title']
    assert len(first.encode('utf-8')) <= 80
    assert body == base_bodies[record['writer']]
    summary['unchanged_non_title_bodies'] += 1
    mol_block = block.split('M  END\n', 1)[0] + 'M  END\n'
    mol = Chem.MolFromMolBlock(mol_block, sanitize=True, removeHs=False, strictParsing=True)
    assert mol is not None
    assert mol.GetProp('_Name') == record['expected_title']
    assert [a.GetAtomicNum() for a in mol.GetAtoms()] == [6, 8]
    assert mol.GetNumBonds() == 1
    assert Chem.MolToSmiles(mol) == 'CO'
    if '> <ORIGINAL_TITLE>' in block:
        supplier = Chem.ForwardSDMolSupplier(io.BytesIO(block.encode('utf-8')), removeHs=False)
        recovered = list(supplier)
        assert len(recovered) == 1 and recovered[0] is not None
        assert recovered[0].GetProp('ORIGINAL_TITLE') == 'retained SD field'
        summary['sdf_property_reads'] += 1
    summary['cases'] += 1
control = subprocess.run(['/opt/homebrew/bin/obabel', '-ismi', '-omol'],
                         input='C ' + '0123456789' * 10 + '\n', text=True,
                         capture_output=True, check=True)
version = subprocess.run(['/opt/homebrew/bin/obabel', '-V'], capture_output=True,
                         text=True, check=True).stdout.strip()
diagnostic = {'version': version, 'input_title_bytes': 100,
              'output_title_bytes': len(control.stdout.splitlines()[0].encode()),
              'stdout': control.stdout, 'stderr': control.stderr,
              'scope': 'Independent synthetic methane title control; not a parity test.'}
(root / 'openbabel-diagnostic.json').write_text(json.dumps(diagnostic, indent=2) + '\n')
(root / 'oracle-summary.json').write_text(json.dumps(summary, indent=2) + '\n')
print(json.dumps(summary, indent=2))
print('Open Babel control title bytes:', diagnostic['output_title_bytes'])
