"""Native PDB fixed-width and CONECT boundary reference inputs."""
import itertools,json
from pathlib import Path
from rdkit import Chem,rdBase
assert rdBase.rdkitVersion=='2026.03.1'

def atom(serial=1,name=' C1 ',element='C',charge='',record='HETATM',residue='UNL',resno=1,xyz=(0.,0.,0.)):
    x,y,z=xyz
    return f'{record:<6}{serial:5d} {name:4} {residue:>3} A{resno:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{1.:6.2f}{0.:6.2f}          {element:>2}{charge:>2}'
blocks=[]
for record,name in itertools.product(['ATOM','HETATM'],[' C1 ','1C  ','CL  ','BR  ','HG  ','HF  ','HO  ','HE  ','HS  ',' Z1 ','    ',' CA ','C   ']):
    blocks.append((f'guess-{record}-{name}',atom(record=record,name=name,element='')+'\nEND\n'))
for length in [6,15,16,17,19,20,21,22,25,26,27,30,37,38,45,46,53,54,59,60,65,66,76,77,78,79,80]:
    blocks.append((f'truncate-{length}',atom(element='Cl',name='CL  ')[:length]+'\nEND\n'))
for charge in ['+1','+2','+ ','++','+0','-1','-2','- ','--','-0',' 1',' -',' +','1 ','2-','0+','x+']:
    blocks.append((f'charge-{charge}',atom(element='N',name=' N  ',charge=charge)+'\nEND\n'))
for off,value in [(30,'-1.25   '),(30,'+.25    '),(30,'12,34   '),(30,'        '),(30,' 1.0e2  '),(54,'  .   '),(60,' +.25 '),(22,' +12'),(6,'12-3 ')]:
    text=atom();text=text[:off]+value+text[off+len(value):]
    blocks.append((f'numeric-{off}-{value}',text+'\nEND\n'))
for endings in ['\n','\r','\r\n','']:
    blocks.append((f'ending-{repr(endings)}',atom(element='C')+endings))
blocks.extend([('endmdl-before-atoms','ENDMDL\n'+atom()+'\nEND\n'),('short-conect',atom()+'\nCONECT\nEND\n')])
base=atom(1,xyz=(0,0,0))+'\n'+atom(2,xyz=(1.5,0,0))+'\n'
for count in [1,2,3,4]:
    for order in [0,1,2]:
        forward='CONECT    1'+'    2'*count+'\n'
        reverse='CONECT    2'+'    1'*count+'\n'
        lines=[forward,reverse] if order==0 else [reverse,forward] if order==1 else [reverse]
        blocks.append((f'conect-{count}-{order}',base+''.join(lines)+'END\n'))
for record in ['CONECTabcde    2','CONECT    1abcde','CONECT    3    1','CONECT    1    3','CONECT    1    1    2','CONECT    1         2']:
    blocks.append((record,base+record+'\nEND\n'))
for residue,resno,element in [('HOH',2,'O'),('UNL',2,'Na'),('ALA',2,'C')]:
    blocks.append((f'cross-residue-{residue}-{element}',atom(1,element='O',name=' O1 ',residue='ALA')+'\n'+atom(2,element=element,residue=residue,resno=resno,xyz=(1.5,0,0))+'\nCONECT    1    2\nEND\n'))
rows=[]
for label,block in blocks:
    mol=Chem.MolFromPDBBlock(block,sanitize=True,removeHs=False,flavor=0,proximityBonding=False)
    unsupported=None
    if mol is not None:
        if any(b.GetBondType()==Chem.BondType.ZERO for b in mol.GetBonds()):unsupported='zero-order bond'
        elif any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()):unsupported='radical'
    row=dict(label=label,block=block,valid=mol is not None,unsupported=unsupported)
    if mol is not None:row.update(atoms=mol.GetNumAtoms(),bonds=mol.GetNumBonds(),smiles=Chem.MolToSmiles(mol))
    rows.append(row)
(Path(__file__).resolve().parents[1]/'validation/rdkit-2026.03.1-pdb-records-boundary.json').write_text(json.dumps(dict(rdkit_version=rdBase.rdkitVersion,rows=rows),indent=2)+'\n')
print(len(rows),'records',sum(r['valid'] and not r['unsupported'] for r in rows),'supported',sum(bool(r['unsupported']) for r in rows),'unsupported')
