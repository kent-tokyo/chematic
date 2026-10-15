#!/usr/bin/env python3
"""Pinned two-way Indigo RDFile interoperability and exact source-field checks."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import chematic
from indigo import Indigo
from rdkit import Chem, rdBase

ROOT=Path(__file__).resolve().parents[1]
FIXTURE=ROOT/'validation/openbabel_issues/rdfile/manifest.json'
INDIGO_VERSION='1.46.0.0-gb36d38f52-aarch64-darwin-appleclang-17000013'


def digest(data):return hashlib.sha256(data).hexdigest()


def mol_key(smiles):
    params=Chem.SmilesParserParams();params.removeHs=False
    m=Chem.MolFromSmiles(smiles,params);assert m is not None
    maps=sorted((a.GetAtomMapNum(),a.GetAtomicNum(),a.GetIsotope(),a.GetFormalCharge()) for a in m.GetAtoms() if a.GetAtomMapNum())
    mapped_bonds=sorted((min(b.GetBeginAtom().GetAtomMapNum(),b.GetEndAtom().GetAtomMapNum()),max(b.GetBeginAtom().GetAtomMapNum(),b.GetEndAtom().GetAtomMapNum()),str(b.GetBondType())) for b in m.GetBonds() if b.GetBeginAtom().GetAtomMapNum() and b.GetEndAtom().GetAtomMapNum())
    return Chem.MolToInchi(m),maps,mapped_bonds


def key(smiles):
    if '>' in smiles:
        parts=smiles.split('>');assert len(parts)==3
        return [sorted(mol_key(s) for s in part.split('.') if s) for part in parts]
    return mol_key(smiles)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--obabel',default='/opt/homebrew/bin/obabel');p.add_argument('--output',type=Path,default=ROOT/'validation/results/openbabel-rdfile-records-2026-10-15.json');args=p.parse_args()
    indigo=Indigo();assert indigo.version()==INDIGO_VERSION,indigo.version()
    assert rdBase.rdkitVersion=='2026.09.1',rdBase.rdkitVersion
    indigo.setOption('molfile-saving-skip-date',True)
    manifest=json.loads(FIXTURE.read_text());payload=FIXTURE.parent/manifest['payload'];text=payload.read_text()
    assert digest(payload.read_bytes())==manifest['sha256']
    base=chematic.parse_rdfile(text);assert base['source_text']==text
    assert [r['kind'] for r in base['records']]==manifest['expected']['kinds']
    projected=chematic.rdfile_to_smiles(text)
    assert len(projected)==5 and projected[-2:]==[None,None]
    assert all(key(a)==key(b) for a,b in zip(projected[:3],manifest['expected']['projection'][:3]))
    assert [f['name'] for f in base['records'][0]['data']]==manifest['expected']['first_field_names']
    assert [f['value'] for f in base['records'][0]['data']]==manifest['expected']['first_field_values']
    values=['','two\n  spaced\n\nlast','long-'+('abc '*70)+'\n  tail\n','internal  whitespace']
    cases=[];forward=0;reverse=0
    with tempfile.TemporaryDirectory(prefix='chematic-rdfile-oracle-') as directory:
        root=Path(directory)
        for order,ending,value in itertools.product(itertools.permutations(range(3)),['\n','\r\n'],values):
            doc={'date_comment':'original synthetic cross-read','records':[]}
            for position,old in enumerate(order):
                r=json.loads(json.dumps(base['records'][old]));r['data']=[{'name':'note','value':value},{'name':'position','value':str(position)}];doc['records'].append(r)
            # Explicit unwrapped profile avoids Indigo's unimplemented '+' joining.
            native=chematic.write_rdfile(doc,wrap_data=False)
            if ending=='\r\n':native=native.replace('\n','\r\n')
            path=root/'native.rdf';path.write_bytes(native.encode())
            raw=chematic.parse_rdfile(native);assert raw['source_text']==native
            assert all(r['data']==d['data'] for r,d in zip(raw['records'],doc['records']))
            objects=list(indigo.iterateRDFile(str(path)));assert len(objects)==3
            native_smiles=chematic.rdfile_to_smiles(native)
            for position,(obj,s,r) in enumerate(zip(objects,native_smiles,doc['records'])):
                assert key(obj.smiles())==key(s)==key(manifest['expected']['projection'][order[position]])
                assert obj.getProperty('note')==value and obj.getProperty('position')==str(position)
                registry='internal-regno' if r['internal_registry'] is not None else 'external-regno'
                expected=r['internal_registry'] or r['external_registry'];assert obj.getProperty(registry)==expected
                forward+=1
            target=root/'indigo.rdf';saver=indigo.writeFile(str(target));saver.rdfHeader()
            for position,old in enumerate(order):
                r=base['records'][old]
                obj=indigo.loadMolecule(r['structure']) if r['kind']=='molecule' else indigo.loadReaction(r['structure'])
                obj.setProperty('note',value);obj.setProperty('position',str(position));saver.rdfAppend(obj)
            saver.close();external=target.read_text();back=chematic.parse_rdfile(external);back_smiles=chematic.rdfile_to_smiles(external)
            assert len(back['records'])==3
            for position,(r,s) in enumerate(zip(back['records'],back_smiles)):
                assert key(s)==key(manifest['expected']['projection'][order[position]])
                props={f['name']:f['value'] for f in r['data']}
                assert props['note']==value and props['position']==str(position)
                reverse+=1
            # Default spec folding remains reversible inside chematic.
            folded=chematic.write_rdfile(doc);restored=chematic.parse_rdfile(folded)
            assert restored['records']==doc['records']
            cases.append({'order':order,'line_ending':'CRLF' if ending=='\r\n' else 'LF','value_sha256':digest(value.encode()),'native_sha256':digest(native.encode()),'indigo_sha256':digest(external.encode()),'folded_sha256':digest(folded.encode())})
        # Record the measured long-data difference separately from passing controls.
        doc={'date_comment':'folding diagnostic','records':[json.loads(json.dumps(base['records'][0]))]}
        long='a'*160+'\n  middle\n';doc['records'][0]['data']=[{'name':'long','value':long}]
        folded=chematic.write_rdfile(doc);path=root/'folded.rdf';path.write_text(folded)
        obj=next(iter(indigo.iterateRDFile(str(path))));observed=obj.getProperty('long')
        assert observed!=long and '+\n' in observed
        diagnostic={'indigo_joins_column_81_plus':False,'expected_value_sha256':digest(long.encode()),'observed_value_sha256':digest(observed.encode()),'folded_file_sha256':digest(folded.encode()),'workaround':'write_rdfile(document, wrap_data=False); explicit plain-long-data extension'}
    version=subprocess.run([args.obabel,'-V'],capture_output=True,text=True,check=True,timeout=10).stdout.strip();assert version.startswith('Open Babel 3.2.1 ')
    formats=subprocess.run([args.obabel,'-L','formats'],capture_output=True,text=True,check=True,timeout=10).stdout
    registered=any(line.split(' -- ')[0].strip() in {'rdf','rdfile'} for line in formats.splitlines())
    assert not registered,'Open Babel RDFile capability changed; inspect before rebaseline'
    failed=subprocess.run([args.obabel,'-irdf','-osmi'],input=text,capture_output=True,text=True,timeout=10)
    assert failed.returncode!=0
    paths=['crates/chematic-mol/src/rdfile.rs','crates/chematic-mol/src/rdfile/reader.rs','crates/chematic-mol/src/rdfile/writer.rs','crates/chematic-py/src/formats/rdfile.rs','crates/chematic-mol/tests/rdfile_contract.rs','crates/chematic-py/tests/test_rdfile.py','crates/chematic-py/python/chematic/__init__.pyi','scripts/validate_rdfile_records.py']
    record={'schema_version':1,'status':'local-source-candidate','issue':'https://github.com/openbabel/openbabel/issues/2681','base':'cac480fb source snapshot plus preserved previous candidates through correlated stereo-group alternatives; patch isolates RDFile records',
        'scope':'RDFile v1 records, V2000/V3000 MOL and V2000 RXN, unresolved registry references and ordered ASCII data; canonical checked writing plus explicit SMILES projection',
        'runtime':{'python':sys.version.split()[0],'indigo':indigo.version(),'rdkit':rdBase.rdkitVersion,'openbabel':version,'chematic_extension_sha256':digest(Path(chematic.__file__).resolve().read_bytes())},
        'source_sha256':{s:digest((ROOT/s).read_bytes()) for s in paths},'fixture':{'path':str(FIXTURE.relative_to(ROOT)),'manifest_sha256':digest(FIXTURE.read_bytes()),**manifest},
        'results':{'datasets':len(cases),'native_to_indigo_records':forward,'indigo_to_native_records':reverse,'total_cross_read_records':forward+reverse,'graph_or_map_mismatches':0,'metadata_mismatches':0,'source_text_mismatches':0,'default_folded_native_roundtrip_datasets':len(cases),'registry_only_records_retained':2,'duplicate_field_entries_retained':2,'case_sha256':digest(json.dumps(cases,sort_keys=True,separators=(',',':')).encode()),'forward_writer_settings':{'wrap_data':False},'reverse_settings':{'molfile-saving-skip-date':True},'graph_comparison':'RDKit InChI by reaction role plus atom-map element/isotope/charge and mapped-bond counters'},
        'indigo_folding_diagnostic':diagnostic,'openbabel_diagnostic':{'rdf_registered':registered,'explicit_read_exit_code':failed.returncode,'stderr':failed.stderr.strip()},
        'limitations':'No embedded V3000 RXN/agents, nested $DATUM $MFMT, both registry identifiers, continued/long headers, non-ASCII RD metadata, query-schema interpretation or registry resolution. Original source/structure payloads remain available. Default output folds ASCII data with column-81 plus; measured Indigo interoperability uses explicit wrap_data=False and supports plain long data. Duplicate metadata is checked against authored ordered vectors, not Indigo properties, which do not preserve duplicates.'}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record['results']))


if __name__=='__main__':main()
