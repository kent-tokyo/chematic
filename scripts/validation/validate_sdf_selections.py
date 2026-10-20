#!/usr/bin/env python3
"""Independently compare indexed SDF selection against COSMolKit 0.5.0."""
import argparse
import hashlib
from importlib.metadata import version
import itertools
import json
from pathlib import Path
import sys
import tempfile

import chematic
import cosmolkit as ck

ROOT=Path(__file__).resolve().parents[2]
FIXTURE=ROOT/'validation/cosmolkit_comparison/sdf_selections/minimal.sdf'
MANIFEST=FIXTURE.with_name('manifest.json')


def digest(data):
    return hashlib.sha256(data).hexdigest()


def keys(n):
    for start,stop,step in itertools.product([None,*range(-8,9)],[None,*range(-8,9)],[-3,-2,-1,1,2,3]):
        key=slice(start,stop,step)
        yield key,list(range(n))[key]
    for count in range(4):
        for key in itertools.product(range(-n,n),repeat=count):
            yield list(key),[i%n for i in key]
    for key in itertools.product([False,True],repeat=n):
        yield list(key),[i for i,b in enumerate(key) if b]
    for key in itertools.permutations(range(n)):
        yield list(key),list(key)


def check_selected(native,indices,expected):
    assert native.source_indices()==indices
    assert native.valid_mask()==[expected['valid_mask'][i] for i in indices]
    assert native.to_smiles_list()==[expected['smiles'][i] for i in indices]
    assert [e['source_index'] for e in native.errors()]==[i for i in indices if not expected['valid_mask'][i]]
    assert [e['position'] for e in native.errors()]==[k for k,i in enumerate(indices) if not expected['valid_mask'][i]]
    for record,index in zip(native.records(),indices):
        if record is not None:
            assert record.name==expected['titles'][index]
            assert record.get('ID')==expected['id_fields'][index]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=ROOT/'validation/results/cosmolkit-sdf-selections-2026-10-15.json')
    args=p.parse_args();assert version('cosmolkit')=='0.5.0'
    payload=FIXTURE.read_bytes();manifest=json.loads(MANIFEST.read_text());assert digest(payload)==manifest['sha256']
    expected=manifest['expected'];n=len(expected['smiles'])
    count=0;successes=0;rejections=0;selected_positions=0;batch_count=0;scalar_count=0
    checksum=hashlib.sha256()
    with tempfile.TemporaryDirectory(prefix='chematic-sdf-oracle-') as tmp:
        for eol in ['LF','CRLF']:
            for terminated in [True,False]:
                text=payload.decode()
                if not terminated:text=text.removesuffix('$$$$\n')
                if eol=='CRLF':text=text.replace('\n','\r\n')
                path=Path(tmp)/'records.sdf';path.write_bytes(text.encode());original=path.read_bytes()
                native=chematic.SdfDataset.open(str(path));oracle=ck.SdfDataset.open(str(path))
                assert len(native)==len(oracle)==n
                for i in range(n):
                    assert native.metadata(i).title==oracle.metadata(i).title()==expected['titles'][i]
                    if expected['valid_mask'][i]:
                        assert native[i].get('ID')==oracle[i].data_field('ID')==expected['id_fields'][i]
                        scalar_count+=1
                for key,indices in keys(n):
                    result=native.select(key,errors='keep');check_selected(result,indices,expected)
                    selected_positions+=len(indices)
                    valid=all(expected['valid_mask'][i] for i in indices)
                    if valid:
                        reference=oracle[key];actual=native[key]
                        assert actual.valid_mask()==reference.valid_mask()
                        assert actual.to_smiles_list()==reference.to_smiles_list()
                        successes+=1
                    else:
                        for operation in [lambda:oracle[key],lambda:native[key]]:
                            try:operation()
                            except ValueError:pass
                            else:raise AssertionError('default indexed selection must raise on malformed input')
                        rejections+=1
                    checksum.update(json.dumps([eol,terminated,repr(key),indices,result.valid_mask(),result.to_smiles_list()],separators=(',',':')).encode())
                    count+=1
                for size in range(1,n+2):
                    batches=list(native.batches(size=size,errors='keep'))
                    reference=list(oracle.batches(size=size,errors='keep'))
                    assert [b.valid_mask() for b in batches]==[b.valid_mask() for b in reference]
                    assert [b.to_smiles_list() for b in batches]==[b.to_smiles_list() for b in reference]
                    assert [i for b in batches for i in b.source_indices()]==list(range(n))
                    batch_count+=len(batches)
                    for b in batches:check_selected(b,b.source_indices(),expected)
                assert path.read_bytes()==original
    paths=['crates/chematic-mol/src/sdf.rs','crates/chematic-mol/src/sdf/selection.rs','crates/chematic-mol/src/sdf/fields.rs','crates/chematic-py/src/io.rs','crates/chematic-py/src/io/dataset.rs','crates/chematic-py/src/io/dataset_selector.rs','crates/chematic-py/src/io/dataset_batch.rs','crates/chematic-py/python/chematic/__init__.pyi','crates/chematic-mol/tests/sdf_selections.rs','crates/chematic-py/tests/test_sdf_selections.py','scripts/validation/validate_sdf_selections.py']
    ck_files=list(Path(ck.__file__).parent.rglob('*.so'))+[Path(ck.__file__)]
    report={'schema_version':1,'status':'local-source-candidate','feature':'indexed SDF selection and position-preserving batch processing',
        'base':'Preserved isolated candidates through GZMAT precision; not integrated into main',
        'runtime':{'python':sys.version.split()[0],'cosmolkit':version('cosmolkit'),'cosmolkit_package_sha256':{str(f.relative_to(Path(ck.__file__).parent)):digest(f.read_bytes()) for f in ck_files},'cosmolkit_settings':'SdfDataset.open defaults; indexing errors raise; batches errors keep','chematic_extension_sha256':digest(Path(chematic.__file__).resolve().read_bytes())},
        'source_sha256':{p:digest((ROOT/p).read_bytes()) for p in paths},'fixture':{'path':str(FIXTURE.relative_to(ROOT)),**manifest},
        'results':{'selection_cases':count,'matched_successful_default_selections':successes,'matched_malformed_default_rejections':rejections,'selected_positions_verified':selected_positions,'batch_values_verified':batch_count,'scalar_metadata_and_field_controls':scalar_count,'variants':4,'selection_order_mismatches':0,'mask_mismatches':0,'smiles_mismatches':0,'field_mismatches':0,'input_mutations':0,'case_sha256':checksum.hexdigest()},
        'scope':'Python slices, integer lists/tuples, strict boolean masks, duplicate/negative indices; keep/raise failure policy; bounded indexed-file batches; mixed V2000/V3000 dispatch',
        'limits':{'default_selection_records':10000,'default_total_source_record_bytes':64<<20,'python_batch_size_max':10000,'count_and_bytes_preflight_before_io':True},
        'limitations':'Uses chematic molecule parsing and native SMILES, not a general COSMolKit/RDKit chemistry parity claim. New SdfSelectionBatch differs from COSMolKit MoleculeBatch. Default indexing raises; select() defaults to keep. Only int/slice/list/tuple selectors; NumPy arrays must be converted to lists. No parallelism/progress UI or generalized MoleculeBatch transforms. Source files must remain unchanged while indexed; I/O/UTF8 failures abort and fuse iteration. Source byte limits bound input work, not exact process RSS. Full BioStructure hierarchy/serialization and broader Layered fingerprint parameter controls remain separate gaps.'}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report['results']))


if __name__=='__main__':main()
