"""Independent CML ID/geometry/schema checks and Open Babel #2120 control."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import xml.etree.ElementTree as ET

import chematic as c
from lxml import etree

p=argparse.ArgumentParser()
p.add_argument('--schema', type=Path, required=True)
p.add_argument('--output', type=Path, required=True)
a=p.parse_args()
root=Path(__file__).resolve().parents[3]
fixtures=Path(__file__).parent
ns={'c':'http://www.xml-cml.org/schema'}
schema=etree.XMLSchema(etree.parse(str(a.schema)))
rows=[]
cases=[('CO',['C_site','O_2']), ('CCO',['C1','C-two','O_third']),
       ('[NH4+]',['N_plus']), ('[13CH3]O',['C_13','O1']),
       ('[2H]O[H]',['D','O','H']), ('C#N',['Carbon','Nitrogen']),
       ('c1ccccc1',[f'Ring_{i}' for i in range(6)])]
for smiles,names in cases:
    mol=c.from_smiles(smiles)
    source_bonds=ET.fromstring(c.write_cml_named_checked(mol,names)).findall('c:bondArray/c:bond',ns)
    expected_bonds=sorted((tuple(sorted(b.get('atomRefs2').split())),b.get('order')) for b in source_bonds)
    permutations=list(itertools.permutations(range(len(names))))[:24]
    for order in permutations:
        permuted=mol.renumber_atoms(list(order))
        ordered_names=[names[i] for i in order]
        for scale in [1.0,1e-200,1e200]:
            coords=[[scale*(i+0.125),-scale*(i+0.25),scale*(i+0.5)] for i in order]
            text=c.write_cml_named_checked(permuted,iter(ordered_names),(row for row in coords))
            schema.assertValid(etree.fromstring(text.encode()))
            xml=ET.fromstring(text)
            atoms=xml.findall('c:atomArray/c:atom',ns)
            assert [x.get('id') for x in atoms]==ordered_names
            assert [[float(x.get(k)) for k in ['x3','y3','z3']] for x in atoms]==coords
            bonds=sorted((tuple(sorted(b.get('atomRefs2').split())),b.get('order')) for b in xml.findall('c:bondArray/c:bond',ns))
            assert bonds==expected_bonds
            parsed=c.from_cml(text)
            # Existing reader checks graph chemistry; it exposes only 2D data.
            assert parsed.formula==mol.formula
            assert parsed.rdkit_smiles==mol.rdkit_smiles
            assert parsed.atom_count==mol.atom_count
            rows.append({'smiles':smiles,'permutation':list(order),'scale':scale,'atom_ids':ordered_names})

source=(fixtures/'water.mol2').read_text()
record=json.loads(c.parse_mol2_record_json(source))
names=[x['name'] for x in record['atoms']]
text=c.write_cml_named_checked(c.from_mol2(source),names,record['coordinates'])
schema.assertValid(etree.fromstring(text.encode()))
water=ET.fromstring(text)
assert [x.get('id') for x in water.findall('c:atomArray/c:atom',ns)]==names
assert [[float(x.get(k)) for k in ['x3','y3','z3']] for x in water.findall('c:atomArray/c:atom',ns)]==record['coordinates']
# Sequential control avoids conflating atom-name replacement with Open Babel's
# handling of nonsequential source IDs.
ob=subprocess.run(['obabel','-imol2',str(fixtures/'water-sequential.mol2'),'-ocml'],capture_output=True,text=True,check=True)
old=ET.fromstring(ob.stdout)
old_ids=[x.get('id') for x in old.findall('c:atomArray/c:atom',ns)]
old_bonds=old.findall('c:bondArray/c:bond',ns)
assert old_ids==['a1','a2','a3'] and len(old_bonds)==2
assert old_ids!=names
# The default chematic writer also generates IDs; this lane is unchanged.
baseline=c.from_mol2(source).to_cml()
assert [x.get('id') for x in ET.fromstring(baseline).findall('c:atomArray/c:atom',ns)]==['a1','a2','a3']
import chematic.chematic as native
report={'source_base':'bef0d008','native_path':native.__file__,
        'native_sha256':hashlib.sha256(Path(native.__file__).read_bytes()).hexdigest(),
        'schema_url':'https://www.xml-cml.org/schema/schema3/schema.xsd',
        'schema_sha256':hashlib.sha256(a.schema.read_bytes()).hexdigest(),
        'independent_schema_geometry_cases':len(rows),'failures':0,
        'mol2_source_ids':[x['source_id'] for x in record['atoms']],
        'mol2_atom_names':names,'candidate_atom_ids':names,
        'openbabel_version':subprocess.check_output(['obabel','-V'],text=True).strip(),
        'openbabel_atom_ids':old_ids,'openbabel_bond_count':len(old_bonds),
        'openbabel_cml':ob.stdout,'openbabel_stderr':ob.stderr,
        'candidate_water_cml':text,'baseline_generated_ids':['a1','a2','a3'],'rows':rows}
a.output.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k not in ['rows','candidate_water_cml','openbabel_cml']},indent=2))
