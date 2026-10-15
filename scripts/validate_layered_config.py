#!/usr/bin/env python3
"""Exact Layered bits/counters against pinned RDKit, with independent packing."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import random
import sys

import chematic
from rdkit import Chem,DataStructs,rdBase,RDLogger

ROOT=Path(__file__).resolve().parents[1]
FIXTURE=ROOT/'validation/rdkit_issues/fingerprints/layered-config.json'


def digest(data):return hashlib.sha256(data).hexdigest()


def pack(bits,size):
    result=bytearray((size+7)//8)
    for bit in bits:result[bit//8]|=1<<(bit%8)
    return bytes(result)


def mask(size,kind,reference):
    if kind=='none':return None
    if kind=='zero':return bytes((size+7)//8)
    if kind=='all':return pack(range(size),size)
    if kind=='alternate':return pack(range(0,size,2),size)
    if kind=='first-hit':return pack(list(reference.GetOnBits())[:1],size)
    if kind=='last':return pack([size-1],size)
    raise AssertionError(kind)


def orders(n):
    if n<2:return [list(range(n))]
    rng=random.Random(9706+n);out=[list(range(n)),list(reversed(range(n)))]
    target=min(6,2 if n==2 else 6)
    while len(out)<target:
        row=list(range(n));rng.shuffle(row)
        if row not in out:out.append(row)
    return out


def settings(n):
    # Every active-layer subset; high unknown flags are ignored by the reference.
    for flags,size in itertools.product([*range(64),0x40,0xffffffff],[1,127,2048]):
        yield {'layer_flags':flags,'nbits':size}
    roots=[None,[]]
    if n:roots += [[0],[n-1],[0,0],[0,n-1]]
    for root,paths,size in itertools.product(roots,[(1,1),(1,2),(2,3),(1,3),(3,7),(1,12)],[9,256,2048]):
        yield {'from_atoms':root,'min_path':paths[0],'max_path':paths[1],'nbits':size,'max_visits':10000000}
        if root:
            yield {'from_atoms':root,'min_path':paths[0],'max_path':paths[1],'nbits':size,'branched_paths':False,'max_visits':10000000}
    # Every supported min/max path combination, including range-dependent closure.
    for lo in range(1,13):
        for hi in range(lo,13):
            yield {'min_path':lo,'max_path':hi,'nbits':257,'max_visits':10000000}
            if n:
                yield {'min_path':lo,'max_path':hi,'nbits':257,'from_atoms':[0],'branched_paths':False,'max_visits':10000000}
    for size in (65535,65536):
        for kind in ('none','all','alternate','first-hit'):
            yield {'nbits':size,'mask_kind':kind,'atom_counts':[23]*n+[123]}
    for flags,size,kind in itertools.product([0,1,4,8,0x25,0x3f],[3,9,255,1024],['zero','all','alternate','first-hit','last']):
        yield {'layer_flags':flags,'nbits':size,'mask_kind':kind,'atom_counts':[2**32-1,7,13][:n]+[3]*max(0,n-3)+[999]}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=ROOT/'validation/results/rdkit-layered-config-2026-10-15.json');args=p.parse_args()
    assert rdBase.rdkitVersion=='2026.09.1';RDLogger.DisableLog('rdApp.warning')
    manifest=json.loads(FIXTURE.read_text());assert digest(json.dumps(manifest['cases'],sort_keys=True,separators=(',',':')).encode())==manifest['sha256']
    checksum=hashlib.sha256();count=0;default_controls=0;counter_values=0;masked=0;linear=0;permuted=0
    for case in manifest['cases']:
        native=chematic.from_smiles(case['smiles']) if case['smiles'] else chematic.RWMol().GetMol();params=Chem.SmilesParserParams();params.removeHs=False
        reference=Chem.MolFromSmiles(case['smiles'],params);assert reference is not None
        n=reference.GetNumAtoms();assert native.atom_count==n
        for order in orders(n):
            m=native.renumber_atoms(order);r=Chem.RenumberAtoms(reference,order) if n else Chem.Mol(reference)
            assert sorted(tuple(sorted(b[:2])) for b in m.bond_table)==sorted(tuple(sorted((b.GetBeginAtomIdx(),b.GetEndAtomIdx()))) for b in r.GetBonds())
            original=m.smiles;m.SetIntProp('trial',9706)
            default=m.rdkit_layered_config();assert default['bits']==bytes(m.rdkit_layered_fp());default_controls+=1
            options=list(settings(n)) if order==list(range(n)) else [
                {'layer_flags':flags,'from_atoms':root,'nbits':size,'branched_paths':branch,'mask_kind':kind,'atom_counts':[17]*n+[123]}
                for flags,root,size,branch,kind in [(0xffffffff,None,2048,True,'none'),(3,[0] if n else None,127,True,'alternate'),(0x25,[n-1] if n else None,9,True,'all'),(8,[0,0] if n else None,256,True,'first-hit')]
            ]
            for original_options in options:
                options=dict(original_options);kind=options.pop('mask_kind','none')
                size=options.get('nbits',2048);seed=options.get('atom_counts',[0]*n);before_seed=list(seed)
                kw={'layerFlags':options.get('layer_flags',0xffffffff),'minPath':options.get('min_path',1),'maxPath':options.get('max_path',7),'fpSize':size,'fromAtoms':options.get('from_atoms'),'branchedPaths':options.get('branched_paths',True)}
                unmasked=Chem.LayeredFingerprint(r,**kw)
                packed_mask=mask(size,kind,unmasked)
                if packed_mask is not None:
                    v=DataStructs.ExplicitBitVect(size)
                    for bit in range(size):
                        if packed_mask[bit//8]&(1<<(bit%8)):v.SetBit(bit)
                    kw['setOnlyBits']=v;options['set_only_bits']=packed_mask;masked+=1
                expected_counts=list(seed);fp=Chem.LayeredFingerprint(r,atomCounts=expected_counts,**kw)
                got=m.rdkit_layered_config(**options)
                expected_bits=pack(fp.GetOnBits(),size)
                assert got['bits']==expected_bits,(case['id'],order,original_options,'bits')
                assert got['atom_counts']==expected_counts,(case['id'],order,original_options,'counts',got['atom_counts'],expected_counts)
                assert seed==before_seed and m.smiles==original and m.GetIntProp('trial')==9706
                assert got['nbits']==size and len(got['bits'])==(size+7)//8
                counter_values+=len(expected_counts);linear+=not kw['branchedPaths'];permuted+=order!=list(range(n))
                checksum.update(json.dumps([case['id'],order,original_options,kind,got['bits'].hex(),got['atom_counts']],sort_keys=True,separators=(',',':')).encode());count+=1
        print(case['id'],count,flush=True)
    paths=['crates/chematic-fp/src/rdkit_layered.rs','crates/chematic-fp/src/rdkit_layered/features.rs','crates/chematic-fp/src/rdkit_layered/config.rs','crates/chematic-fp/src/rdkit_layered/config/validation.rs','crates/chematic-fp/src/rdkit_layered/config/enumeration.rs','crates/chematic-fp/src/rdkit_layered/config/enumeration/linear.rs','crates/chematic-py/src/mol_methods/layered_config.rs','crates/chematic-py/python/chematic/__init__.pyi','crates/chematic-fp/tests/rdkit_layered_config.rs','crates/chematic-py/tests/test_layered_config.py','scripts/validate_layered_config.py']
    report={'schema_version':1,'status':'local-source-candidate','reference':'RDKit LayeredFingerprint 2026.09.1; experimental family',
        'base':'Preserved isolated candidates through indexed SDF selections; no main integration',
        'runtime':{'python':sys.version.split()[0],'rdkit':rdBase.rdkitVersion,'settings':'removeHs=False; use indexed explicit-H graphs; all masks constructed as ExplicitBitVect; original molecule properties retained','chematic_extension_sha256':digest(Path(chematic.__file__).resolve().read_bytes())},
        'source_sha256':{p:digest((ROOT/p).read_bytes()) for p in paths},'fixture':{'path':str(FIXTURE.relative_to(ROOT)),'manifest_sha256':digest(FIXTURE.read_bytes()),**manifest},
        'results':{'cases':count,'fixtures':len(manifest['cases']),'all_active_layer_subsets':64,'all_supported_path_ranges':78,'default_legacy_controls':default_controls,'seeded_counter_values_verified':counter_values,'masked_cases':masked,'rooted_linear_cases':linear,'permuted_cases':permuted,'bits_mismatches':0,'counter_mismatches':0,'input_mutations':0,'case_sha256':checksum.hexdigest()},
        'limits':{'atoms':512,'bonds':1024,'maximum_path_bonds':12,'maximum_bits':65536,'default_visits':1000000,'maximum_visits':10000000,'default_selected_paths':100000,'maximum_selected_paths':1000000,'linear_frontier_uses_path_cap':True},
        'limitations':'Checked concrete molecules only: no complex query bonds, dative/zero bonds, loops/parallel edges or index-changing sanitized views. Unrooted branchedPaths=False is explicitly unsupported: the pinned RDKit source passes atom paths as bond indices and isolated CC/CCO probes in this session exited with signal 11. Rooted linear preserves the source range-wide final-step closure convention. No general chemical parity, unrestricted parameter range, performance or release claim. Seed counters are copied, not mutated in place; u32 wrapping and trailing entries are preserved. Search visits and selected_paths are implementation diagnostics, not RDKit API fields.'}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report['results']))


if __name__=='__main__':main()
