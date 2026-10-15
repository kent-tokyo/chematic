#!/usr/bin/env python3
"""Independent GRO field checks and pinned Open Babel #2817 diagnostic."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import chematic

ROOT=Path(__file__).resolve().parents[1]
FIXTURE=ROOT/'validation/openbabel_issues/gro/two-letter-elements.gro'
SPEC='https://manual.gromacs.org/current/reference-manual/file-formats.html#gro'


def digest(data):
    return hashlib.sha256(data).hexdigest()


def xyz(text):
    lines=text.splitlines();count=int(lines[0])
    assert len(lines)>=count+2
    return [(v[0],[float(n) for n in v[1:4]]) for v in (s.split() for s in lines[2:count+2])]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--obabel',default='/opt/homebrew/bin/obabel')
    parser.add_argument('--output',type=Path,default=ROOT/'validation/results/openbabel-gro-elements-2026-10-15.json')
    args=parser.parse_args()
    version=subprocess.run([args.obabel,'-V'],capture_output=True,text=True,check=True,timeout=10).stdout.strip()
    if not version.startswith('Open Babel 3.2.1 '):raise RuntimeError(f'expected Open Babel 3.2.1, got {version}')
    text=FIXTURE.read_text()
    unresolved=chematic.parse_gro(text)
    resolved=chematic.parse_gro(text,elements=['Be','Ce','Fe'])
    assert [a['element'] for a in unresolved['atoms']]==[None]*3
    assert [a['element'] for a in resolved['atoms']]==['Be','Ce','Fe']
    process=subprocess.run([args.obabel,'-igro',str(FIXTURE),'-oxyz'],capture_output=True,text=True,check=True,timeout=10)
    output=xyz(process.stdout)
    assert [a[0] for a in output]==['B','C','F'], 'upstream behavior changed; inspect before rebaseline'
    assert [a[1] for a in output]==[[0,0,0],[20,0,0],[40,0,0]]
    # Explicit element-column control, not a new GRO-to-XYZ product adapter.
    control='3\noriginal element assignments, positions scaled nm to angstrom\n'+'\n'.join(a['element']+' '+' '.join(str(v*10) for v in a['position_nm']) for a in resolved['atoms'])+'\n'
    with tempfile.TemporaryDirectory(prefix='chematic-gro-control-') as tmp:
        path=Path(tmp)/'explicit.xyz';path.write_text(control)
        cross=subprocess.run([args.obabel,'-ixyz',str(path),'-oxyz'],capture_output=True,text=True,check=True,timeout=10)
    assert [a[0] for a in xyz(cross.stdout)]==['Be','Ce','Fe']
    assert [a[1] for a in xyz(cross.stdout)]==[[0,0,0],[20,0,0],[40,0,0]]
    names=['BE','CE','FE','CA','CL','NA'];elements=['Be','Ce','Fe','C','Cl','Na']
    cases=[]
    for case,order in enumerate(itertools.permutations(range(6))):
        precision=(3,5,8)[case%3];width=precision+5
        velocity=bool((case//3)%2);triclinic=bool((case//6)%2)
        lines=['original synthetic six-name profile','6']
        for old in order:
            line=f'{old+1:5d}{"UNL":<5}{names[old]:>5}{old+1:5d}'
            line+=''.join(f'{v:{width}.{precision}f}' for v in [old*2+.125,-.25,.5])
            if velocity:line+=''.join(f'{v:{width}.{precision+1}f}' for v in [-.125,.25,.5])
            lines.append(line)
        lines.append('8 9 10 0 0 .4 0 .5 .6' if triclinic else '8 9 10')
        source='\n'.join(lines)+'\n'
        raw=chematic.parse_gro(source)
        d=chematic.parse_gro(source,elements=[elements[i] for i in order])
        assert raw['source_text']==d['source_text']==source
        assert raw['element_policy']=='unresolved' and d['element_policy']=='caller_supplied'
        assert all(a['element'] is None for a in raw['atoms'])
        for i,old in enumerate(order):
            a=d['atoms'][i]
            # Known generator width is an independent fixed-column oracle;
            # it does not reuse the reader's decimal-distance inference.
            source_line=lines[i+2]
            values=[float(source_line[20+j*width:20+(j+1)*width]) for j in range(3)]
            assert a['position_nm']==values==[old*2+.125,-.25,.5]
            assert a['atom_name']==source_line[10:15].strip()==names[old]
            assert a['residue_number']==int(source_line[:5])==old+1
            assert a['atom_number']==int(source_line[15:20])==old+1
            assert a['residue_name']=='UNL' and a['element']==elements[old]
            assert a['velocity_nm_per_ps']==([-.125,.25,.5] if velocity else None)
        expected_box=[[8,0,0],[.4,9,0],[.5,.6,10]] if triclinic else [[8,0,0],[0,9,0],[0,0,10]]
        assert d['box_vectors_nm']==expected_box
        assert len(d['box_values_nm'])==(9 if triclinic else 3)
        cases.append({'order':order,'precision':precision,'velocity':velocity,'triclinic':triclinic,'source_sha256':digest(source.encode())})
    paths=['crates/chematic-mol/src/gro.rs','crates/chematic-py/src/formats/gro.rs','crates/chematic-mol/tests/gro_contract.rs','crates/chematic-py/tests/test_gro.py','crates/chematic-py/python/chematic/__init__.pyi','scripts/validate_gro_elements.py']
    record={
        'schema_version':1,'status':'local-source-candidate','issue':'https://github.com/openbabel/openbabel/issues/2817',
        'base':'cac480fb source snapshot plus preserved earlier RDKit/Open Babel candidates; patch isolates GRO reader',
        'scope':'single-frame GRO record reader and explicit element resolution; no guessed elements, bonds, topology or new graph converter',
        'runtime':{'python':sys.version.split()[0],'openbabel':version,'chematic_extension_sha256':digest(Path(chematic.__file__).resolve().read_bytes())},
        'source_sha256':{s:digest((ROOT/s).read_bytes()) for s in paths},
        'fixture':{'id':'synthetic-be-ce-fe-gro-v1','path':str(FIXTURE.relative_to(ROOT)),'sha256':digest(text.encode()),
            'provenance':{'source':'original synthetic fixed-column rows authored from the GROMACS GRO field definition','date':'2026-10-15','spec_retrieved':'2026-10-15','spec_url':SPEC,'external_data':'no upstream attachment, coordinates or POSCAR copied','transformations':'three separated atoms named BE/CE/FE; explicit Be/Ce/Fe caller assignment; generated six-name validation adds ambiguous CA, CL and NA with source-authored C, Cl and Na assignments'},
            'license':'MIT OR Apache-2.0','redistribution':{'status':'allowed','rationale':'original project fixture'},
            'minimization':'three rows separate first-letter truncation from chemistry; no bonds or molecule reconstruction',
            'oracle':{'name':'source-authored elements and independent fixed-column decoder','version':'GROMACS 2026.4 GRO manual; Open Babel 3.2.1 defect/control diagnostic','url':SPEC,'settings':'nm and nm/ps retained; variable precision and 3/9-value box profiles; no atom-name heuristic'},
            'expected':{'unresolved_elements':[None,None,None],'explicit_elements':['Be','Ce','Fe'],'coordinate_unit':'nm','positions_nm':[[0,0,0],[2,0,0],[4,0,0]],'numeric_tolerance':0.0,'fixture_encoding':'UTF-8 LF'}},
        'results':{'permutations':len(cases),'field_mismatches':0,'element_mismatches':0,'source_text_mismatches':0,'case_sha256':digest(json.dumps(cases,sort_keys=True,separators=(',',':')).encode()),'precisions':[3,5,8],'velocity_presence':[False,True],'box_component_counts':[3,9]},
        'openbabel':{'command':['obabel','-igro','<fixture>','-oxyz'],'output_elements':[a[0] for a in output],'output_sha256':digest(process.stdout.encode()),'stderr':process.stderr.strip(),'explicit_element_XYZ_control':{'input_sha256':digest(control.encode()),'output_sha256':digest(cross.stdout.encode()),'output_elements':[a[0] for a in xyz(cross.stdout)],'passed':True}},
        'limits':'One frame only; fixed five-character identifiers, n+5 numeric fields with n>=1, optional same-width n+1-decimal velocities, finite values and default input/atom/line limits. Raw box vectors are preserved without periodicity, invertibility or simulation-cell validation. No topology, charges, isotopes, writer or trajectory API.',
    }
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record['results']))


if __name__=='__main__':main()
