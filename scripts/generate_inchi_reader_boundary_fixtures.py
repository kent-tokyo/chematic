"""Build references from RDKit's fixed-version official InChI regressions.

Pass External/INCHI-API/test.cpp from the Release_2026_03_1 tag. Its hash
is checked before extracting adjacent C++ string literals; no implementation
output is used to select cases or expected values.
"""
import ast
import hashlib
from pathlib import Path
import re
import sys
from rdkit import Chem, rdBase

assert rdBase.rdkitVersion == '2026.03.1', rdBase.rdkitVersion
source = Path(sys.argv[1]).read_bytes()
assert hashlib.sha256(source).hexdigest() == 'd3faba4a81f82c03647c275f0c0978cd7887e883b0fd9a3ca3e473252d1a9b01'
literal = r'"(?:[^"\\]|\\.)*"'
values = []
for match in re.finditer(literal + r'(?:\s*' + literal + r')*', source.decode()):
    if not match.group().startswith('"InChI='):
        continue
    value = ''.join(ast.literal_eval(part) for part in re.findall(literal, match.group()))
    if value.startswith('InChI=') and value not in values:
        values.append(value)
assert len(values) == 37
rows = []
for inchi in values:
    mol = Chem.MolFromInchi(inchi)
    expected = Chem.MolToSmiles(mol) if mol is not None else '<NONE>'
    rows.append(inchi + '\t' + expected)
root = Path(__file__).resolve().parents[1]
(root / 'validation/rdkit-2026.03.1-inchi-reader-boundary.tsv').write_text('\n'.join(rows) + '\n')
print(f'wrote {len(rows)} cases')
