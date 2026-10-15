#!/usr/bin/env python3
"""Independently verify graph witnesses and allowed configurations for group reduction."""
import argparse
import hashlib
import json
from pathlib import Path
import random
import sys
import chematic
from rdkit import Chem,RDLogger,rdBase
from rdkit.Chem import rdDepictor

ROOT=Path(__file__).resolve().parents[1]
FIXTURE=ROOT/'validation/rdkit_issues/stereo/meso-group-proof.json'


def digest(data):return hashlib.sha256(data).hexdigest()


def parse(s):
    p=Chem.SmilesParserParams();p.removeHs=False;m=Chem.MolFromSmiles(s,p);assert m is not None;return m


def inchi(s):return Chem.MolToInchi(parse(s))


def reference(base,native,specs,mask):
    m=Chem.Mol(base)
    for bit,(_,_,atoms,bonds) in enumerate(g for g in specs if g[0]!='abs'):
        if not mask&(1<<bit):continue
        for idx in atoms:m.GetAtomWithIdx(idx).InvertChirality()
        for idx in bonds:
            a,b,*_=native.bond_table[idx];bond=m.GetBondBetweenAtoms(a,b);stereo=bond.GetStereo()
            assert stereo in (Chem.BondStereo.STEREOE,Chem.BondStereo.STEREOZ)
            bond.SetStereo(Chem.BondStereo.STEREOZ if stereo==Chem.BondStereo.STEREOE else Chem.BondStereo.STEREOE)
    for b in m.GetBonds():b.SetBondDir(Chem.BondDir.NONE)
    Chem.SetDoubleBondNeighborDirections(m)
    assert inchi(Chem.MolToSmiles(m))==Chem.MolToInchi(m)
    return m


def orders(n):
    rng=random.Random(9705);out=[list(range(n)),list(reversed(range(n)))]
    while len(out)<24:
        order=list(range(n));rng.shuffle(order)
        if order not in out:out.append(order)
    return out


def graph_signature(m):
    return [(a.GetAtomicNum(),a.GetIsotope(),a.GetFormalCharge(),a.GetIsAromatic(),a.GetAtomMapNum()) for a in m.GetAtoms()]


def verify_witness(a,b,mapping,native):
    assert sorted(mapping)==list(range(a.GetNumAtoms()))
    matches=b.GetSubstructMatches(a,uniquify=False,useChirality=True,maxMatches=200000)
    assert tuple(mapping) in matches,'native witness is not an RDKit chirality-aware complete graph match'
    assert [graph_signature(b)[j] for j in mapping]==graph_signature(a)
    for _,_,atoms,bonds in native.stereo_groups:
        assert sorted(mapping[i] for i in atoms)==sorted(atoms),'witness exchanges enhanced-group membership'
        endpoints={tuple(sorted(native.bond_table[i][:2])) for i in bonds}
        mapped={tuple(sorted((mapping[x],mapping[y]))) for x,y in endpoints}
        assert mapped==endpoints,'witness exchanges independent bond groups'


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'validation/results/rdkit-meso-stereo-groups-2026-10-15.json');args=p.parse_args()
    assert rdBase.rdkitVersion=='2026.09.1';Chem.SetUseLegacyStereoPerception(True);RDLogger.DisableLog('rdApp.warning')
    fixture=json.loads(FIXTURE.read_text());assert digest(json.dumps(fixture['cases'],sort_keys=True,separators=(',',':')).encode())==fixture['sha256']
    checked=[];witness_count=0;positive=0;negative=0;roundtrips=0;contexts=0
    for case in fixture['cases']:
        base=chematic.from_smiles(case['smiles']);rd=parse(case['smiles']);assert base.atom_count==rd.GetNumAtoms()
        assert [(a[1],a[2],a[3]) for a in base.atom_table]==[(a.GetAtomicNum(),a.GetFormalCharge(),a.GetIsAromatic()) for a in rd.GetAtoms()]
        norm=lambda bs:sorted((min(a,b),max(a,b),str(kind),aromatic) for a,b,kind,aromatic in bs)
        assert norm(base.bond_table)==norm([(b.GetBeginAtomIdx(),b.GetEndAtomIdx(),b.GetBondType(),b.GetIsAromatic()) for b in rd.GetBonds()])
        drawing=Chem.Mol(rd);rdDepictor.Compute2DCoords(drawing);coords=drawing.GetConformer().GetPositions()[:,:2].tolist()
        group_orders=[list(range(len(case['groups'])))]
        if len(case['groups'])>1:group_orders.append(list(reversed(group_orders[0])))
        for kind in ('and','or'):
            for group_order in group_orders:
                specs=[list(case['groups'][i]) for i in group_order];selected=group_order.index(case['selected']);specs[selected][0]=kind
                relative=[i for i,g in enumerate(specs) if g[0]!='abs'];bit=relative.index(selected)
                refs=[reference(rd,base,specs,mask) for mask in range(1<<len(relative))]
                allowed={Chem.MolToInchi(m) for m in refs}
                for order in orders(base.atom_count):
                    source=base.with_stereo_groups([tuple(g) for g in specs]).renumber_atoms(order);before=source.stereo_groups;original=source.smiles
                    source.SetIntProp('trial',9705);report=source.reduce_meso_stereo_group(selected);out=report['molecule']
                    assert report['reduced']==case['expected_reduced']
                    assert out.atom_table==source.atom_table and out.bond_table==source.bond_table and inchi(out.smiles)==inchi(original)
                    assert out.GetIntProp('trial')==9705 and source.stereo_groups==before and source.smiles==original
                    _,family=out.enumerate_stereo_group_alternatives();assert {inchi(m.smiles) for m in family}==allowed
                    if report['reduced']:
                        positive+=1;assert len(report['witnesses'])==len(refs)//2==report['checked_pairs']
                        expected=[];merged=False
                        selected_atoms=before[selected][2];selected_bonds=before[selected][3]
                        for i,(k,num,atoms,bonds) in enumerate(before):
                            if i==selected:continue
                            if k=='abs':atoms=atoms+selected_atoms;bonds=bonds+selected_bonds;merged=True
                            expected.append((k,num,atoms,bonds))
                        if not merged:
                            expected=list(before);expected[selected]=('abs',0,selected_atoms,selected_bonds)
                        assert out.stereo_groups==expected
                        for mask,mapping in report['witnesses']:
                            a=Chem.RenumberAtoms(refs[mask],order);b=Chem.RenumberAtoms(refs[mask|(1<<bit)],order)
                            verify_witness(a,b,mapping,source);witness_count+=1
                    else:
                        negative+=1;assert out.stereo_groups==before
                    block=out.to_mol_v3000([coords[i] for i in order],strict=True);back=chematic.from_mol_v3000(block)
                    assert back.stereo_groups==out.stereo_groups
                    assert inchi(back.smiles)==inchi(out.smiles)
                    roundtrips+=1;contexts+=report['checked_pairs']
                    checked.append({'id':case['id'],'kind':kind,'group_order':group_order,'order':order,'reduced':report['reduced'],'checked_pairs':report['checked_pairs'],'visits':report['search_visits'],'witnesses':report['witnesses']})
    # Separate upstream canonicalization diagnostics and atom-only control.
    def upstream(s,atoms,bonds):
        rw=Chem.RWMol(parse(s));rw.SetStereoGroups([Chem.CreateStereoGroup(Chem.StereoGroupType.STEREO_AND,rw,atoms,bonds)])
        m=Chem.CanonicalizeStereoGroups(rw.GetMol(),Chem.StereoGroupAbsOptions.AlwaysInclude)
        return [{'kind':g.GetGroupType().name,'atoms':[a.GetIdx() for a in g.GetAtoms()],'bonds':[b.GetIdx() for b in g.GetBonds()]} for g in m.GetStereoGroups()]
    mixed=upstream(r'C/C=C/[C@H](C)/C=C\C',[3],[1,5]);assert sum(len(g['bonds']) for g in mixed)==0
    atom_control=upstream('C[C@H](O)[C@H](O)C',[1,3],[]);assert all(g['kind']=='STEREO_ABSOLUTE' for g in atom_control)
    paths=['crates/chematic-chem/src/meso_stereo_group.rs','crates/chematic-chem/src/meso_stereo_group/graph.rs','crates/chematic-chem/src/meso_stereo_group/search.rs','crates/chematic-py/src/mol_methods/stereo_groups.rs','crates/chematic-chem/tests/meso_stereo_group.rs','crates/chematic-py/tests/test_meso_stereo_group.py','crates/chematic-py/python/chematic/__init__.pyi','scripts/validate_meso_stereo_groups.py']
    result={'schema_version':1,'status':'local-source-candidate','issue':'https://github.com/rdkit/rdkit/issues/9705','base':'cac480fb source snapshot plus preserved previous candidates through RDFile; requires checked groups, full inversion and correlated alternatives',
        'scope':'conservative proof-backed promotion of one AND/OR group to ABS; exact group/label-colored formal graph mappings in every other-group assignment, no general canonicalization',
        'runtime':{'python':sys.version.split()[0],'rdkit':rdBase.rdkitVersion,'settings':'legacy stereo; retain explicit H; direct atom/bond flips; chirality-aware complete graph mappings with separate group/label checks','chematic_extension_sha256':digest(Path(chematic.__file__).resolve().read_bytes())},
        'source_sha256':{s:digest((ROOT/s).read_bytes()) for s in paths},'fixture':{'path':str(FIXTURE.relative_to(ROOT)),'manifest_sha256':digest(FIXTURE.read_bytes()),**fixture},
        'results':{'cases':len(checked),'fixtures':len(fixture['cases']),'promoted':positive,'retained':negative,'witnesses_verified_by_rdkit':witness_count,'configuration_context_pairs':contexts,'v3000_roundtrips':roundtrips,'decision_mismatches':0,'witness_mismatches':0,'allowed_configuration_set_mismatches':0,'metadata_mismatches':0,'input_mutations':0,'numberings_per_fixture_kind_group_order':24,'case_sha256':digest(json.dumps(checked,sort_keys=True,separators=(',',':')).encode())},
        'rdkit_mixed_canonicalization_diagnostic':{'input_atoms':[3],'input_bonds':[1,5],'output_groups':mixed,'bond_members_dropped':True},'rdkit_atom_only_meso_control':{'output_groups':atom_control,'absolute_only':True},
        'limits':{'atoms':64,'bonds':128,'relative_groups':6,'family_representatives':64,'default_candidate_visits':1000000,'maximum_candidate_visits':10000000,'budget_unknown_is_error':True},
        'limitations':'Preserving group identities, atom maps and tags can conservatively retain a chemically redundant group. Retained does not mean a general non-meso classification. Undefined locally tetrahedral sites, wildcard/R groups, non-tetrahedral/allene/axial stereo, query/dative/zero bonds, loops/parallel edges, stashed/conflicting/unaccounted directions fail closed. No 3D feasibility, stereogenicity cleanup, canonical representative, group renumbering, CXSMILES extension or general molecular identity claim. Same-number mixed-group interpretation follows #9705 proposal, not confirmed BIOVIA interoperability.'}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result['results']))


if __name__=='__main__':main()
