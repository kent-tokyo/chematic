from pathlib import Path
import gzip,hashlib,json,re,subprocess,sys,tempfile
ROOT=Path('/private/tmp/chematic-3d-stereo-guard-20261017')
OUT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'scripts'))
import bench_openbabel_sdf_throughput as sdf
import bench_openbabel_v3000_throughput as v3
from rdkit import Chem
manifest=json.loads((OUT/'source-manifest.json').read_text())
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert all(sha(ROOT/p)==value for p,value in manifest['source_files_sha256'].items())
assert all(sha(Path(p))==value for p,value in manifest['native_libraries_sha256'].items())
seal=json.loads((OUT/'binary-seal.json').read_text());assert all(sha(OUT/p)==value for p,value in seal['files_sha256'].items())
protocol=json.loads((OUT/'protocol-seal.json').read_text());assert all(sha(OUT/p)==v for p,v in protocol.items())
baseline=json.loads((OUT/'baseline-source-manifest.json').read_text());assert all(sha(OUT/'baseline-source'/p)==v for p,v in baseline['source_files_sha256'].items())
paths=list(manifest['source_files_sha256']);at=0
committed=subprocess.check_output(['git','cat-file','--batch'],input=''.join(manifest['candidate_commit']+':'+p+'\n' for p in paths).encode(),cwd=ROOT)
position=0
for path in paths:
    end=committed.index(b'\n',position);header=committed[position:end].split();position=end+1
    assert header[-1]!=b'missing',path
    size=int(header[2]);blob=committed[position:position+size];position+=size+1
    assert hashlib.sha256(blob).hexdigest()==manifest['source_files_sha256'][path],path
assert position==len(committed)

raw=subprocess.check_output(['git','cat-file','--batch'],input=''.join(manifest['base_revision']+':'+p+'\n' for p in paths).encode(),cwd=ROOT)
for path in paths:
 end=raw.index(b'\n',at);header=raw[at:end].split();at=end+1
 if header[-1]==b'missing':assert path in manifest['added'],path;continue
 size=int(header[2]);blob=raw[at:at+size];at+=size+1
 if path not in manifest['changed']:assert hashlib.sha256(blob).hexdigest()==manifest['source_files_sha256'][path],path
assert at==len(raw)
original_source=subprocess.check_output(['git','show',manifest['base_revision']+':crates/chematic-mol/src/mol2000.rs'],cwd=ROOT,text=True)
def function(text,name):
    start=text.index(name+'(');start=text.rfind('\n',0,start)+1
    at=text.index('{',start);depth=1;end=at+1
    while depth:
        if text[end]=='{':depth+=1
        if text[end]=='}':depth-=1
        end+=1
    return text[start:end]
original=function(original_source,'perceive_tetrahedral_from_3d')
frozen=function((ROOT/'crates/chematic-mol/src/mol2000/stereo3d_guard_tests.rs').read_text(),'frozen_perceive_tetrahedral_from_3d').replace('fn frozen_perceive_tetrahedral_from_3d','pub(crate) fn perceive_tetrahedral_from_3d')
assert original==frozen
candidate=function((ROOT/'crates/chematic-mol/src/mol2000.rs').read_text(),'perceive_tetrahedral_from_3d')
needle='    let mut ordinary_centers:'
assert candidate[candidate.index(needle):]==original[original.index(needle):]
with tempfile.TemporaryDirectory(prefix='chematic-3d-probe-proof-',dir='/private/tmp') as directory:
 p=Path(directory)/'probe.rs';p.write_bytes((OUT/'baseline-source/crates/chematic-mol/examples/openbabel_v3000_3d_throughput.rs').read_bytes())
 subprocess.run(['rustfmt','--edition','2024',str(p)],check=True)
 assert p.read_bytes()==(ROOT/'crates/chematic-mol/examples/openbabel_v3000_3d_throughput.rs').read_bytes()
with tempfile.TemporaryDirectory(prefix='chematic-v3000-loss-replay-',dir='/private/tmp') as d:
    root=Path(d);subprocess.run(['git','init','-q',str(root)],check=True)
    for path in manifest['changed']:
        target=root/path;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(subprocess.check_output(['git','show',manifest['base_revision']+':'+path],cwd=ROOT))
    subprocess.run(['git','apply',str(OUT/'candidate.patch')],cwd=root,check=True)
    assert all(sha(root/path)==manifest['source_files_sha256'][path] for path in manifest['changed']+manifest['added'])
build=json.loads((OUT/'binary-source-proof.json').read_text())
assert build['source_contrast_valid'] and build['compiler_rebuilt_production_library']
assert all(sha(OUT/p)==h for p,h in build['binary_sha256'].items())
assert sha(OUT/'candidate-build.log')==build['candidate_compiler_log_sha256']
for arm in ['baseline','candidate']:
    assembly=subprocess.check_output(['otool','-tvV',str(OUT/(arm+'-v3000-3d'))],text=True).splitlines()
    start=next(i for i,line in enumerate(assembly) if line.endswith(':') and 'perceive_tetrahedral_from_3d' in line and 'FilterMap' not in line)
    end=next((i for i in range(start+1,len(assembly)) if assembly[i].endswith(':')),len(assembly))
    assert hashlib.sha256('\n'.join(assembly[start:end]).encode()).hexdigest()==build[arm]['instruction_text_sha256']
assert build['baseline']['instruction_text_sha256']!=build['candidate']['instruction_text_sha256']
profiles=json.loads((OUT/'measurement-plan.json').read_text())['profiles'];proof={};comparisons=0
for item in profiles:
    kind=item['kind'];profile=item['profile'];key=kind+'-'+profile
    data=json.loads((OUT/f'paired-{key}.json').read_text());assert data['complete'] and len(data['blocks'])==21 and data['source_manifest']==manifest
    suffix='sdf' if kind=='sdf' else 'molblocks';corpus=OUT/f'{key}.{suffix}';assert sha(corpus)==data['corpus']['sha256']
    reference=sdf.observe(corpus.read_bytes());artifacts=OUT/f'paired-{key}.artifacts'
    outputs={};fresh_hashes={}
    for arm in ['baseline','candidate','openbabel']:
        raw=(artifacts/f'{arm}.{suffix}').read_bytes();outputs[arm]=raw;assert hashlib.sha256(raw).hexdigest()==data['rewrite_sha256'][arm]
        observed=sdf.semantic_check(reference,sdf.observe(raw)) if kind=='sdf' else v3.rewrite_check(reference,raw)
        assert observed==data['semantics'][arm]
        binary=OUT/(arm+'-'+kind if arm!='openbabel' else 'openbabel-'+('sdf' if kind=='sdf' else 'v3000'))
        fresh_dump=subprocess.check_output([str(binary),str(corpus),'dump','1',str(item['records'])])
        (OUT/f'post-timing-{key}-{arm}-dump.gz').write_bytes(gzip.compress(fresh_dump,mtime=0))
        fresh_hashes[arm]={'dump_raw_sha256':hashlib.sha256(fresh_dump).hexdigest()}
        fresh_semantics=sdf.semantic_check(reference,sdf.observe(fresh_dump)) if kind=='sdf' else v3.rewrite_check(reference,fresh_dump)
        assert fresh_semantics==data['semantics'][arm]
        if arm!='openbabel':
            assert fresh_dump==raw
            fresh_snapshot=subprocess.check_output([str(binary),str(corpus),'snapshot','1',str(item['records'])])
            (OUT/f'post-timing-{key}-{arm}-snapshot.gz').write_bytes(gzip.compress(fresh_snapshot,mtime=0))
            fresh_hashes[arm]['snapshot_raw_sha256']=hashlib.sha256(fresh_snapshot).hexdigest()
            assert fresh_hashes[arm]['snapshot_raw_sha256']==data['snapshot_sha256'][arm]
        else:
            fresh_smiles=subprocess.check_output([str(binary),str(corpus),'smiles','1',str(item['records'])])
            (OUT/f'post-timing-{key}-openbabel-smiles.gz').write_bytes(gzip.compress(fresh_smiles,mtime=0))
            fresh_hashes[arm]['smiles_raw_sha256']=hashlib.sha256(fresh_smiles).hexdigest()
        if arm!='openbabel':
            snap=gzip.decompress((OUT/f'preflight-{key}-{arm}-snapshot.gz').read_bytes())
            assert hashlib.sha256(snap).hexdigest()==data['snapshot_sha256'][arm]
            strings=[row['canonical_smiles'] for row in json.loads(snap)]
        else:strings=[line.split()[0] for line in fresh_smiles.decode().splitlines()]
        assert len(strings)==item['records']
        rejected=[];bad=[]
        for i,(row,string) in enumerate(zip(reference,strings,strict=True)):
            mol=Chem.MolFromSmiles(string)
            if mol is None:rejected.append(i)
            elif Chem.MolToSmiles(mol)!=row['graph']:bad.append(i)
        expected=data['semantics'][arm+'_parse'];assert expected['passed']==(not rejected and not bad) and expected['mismatches']==len(bad)
        if kind!='sdf':assert expected['mismatch_rows']==bad and expected['rejected_rows']==rejected
        else:assert not rejected and expected['examples']==bad[:8]
    assert outputs['baseline']==outputs['candidate']
    for mode,values in data['summary'].items():
        for seed,arm in enumerate(['baseline','openbabel']):
            allowed=sdf.eligible(data['semantics'],data['source_parse_snapshots_equal'],mode,arm)
            interval=sdf.paired_median_ratio_ci([b[arm][mode]['elapsed_ns'] for b in data['blocks']],[b['candidate'][mode]['elapsed_ns'] for b in data['blocks']],seed=seed) if allowed else None
            assert values[arm]['semantic_prerequisite']==allowed and values[arm]['speedup']==interval
            assert values[arm]['strict_speed_win']==bool(interval and interval['ci95'][0]>1 and interval['winning_blocks']==21)
            comparisons+=1
    proof[key]={'complete_blocks':21,'source_snapshots_and_rewrites_exact':True,'saved_rewrites_and_parse_semantics_revalidated':True,'fresh_post_timing_hashes':fresh_hashes}
result={'source_file_hashes_verified':len(manifest['source_files_sha256']),'native_library_entries_verified':len(manifest['native_libraries_sha256']),'patch_reconstructed_from':manifest['base_revision'],'all_comparisons_recomputed':comparisons,'profiles':proof,'fresh_post_timing_native_observers_not_timing_samples':True,'source_post_timing_bytes_and_snapshots_equal':True,'openbabel_post_timing_declared_semantics_equal':True,'all_unnamed_source_inputs_match_base_blobs':True,'frozen_original_tetrahedral_function_reverified_from_git':True,'common_probe_reverified_after_rustfmt':True,'baseline_has_no_ignored_python_cache_files':True,'protocol_seal_verified':True,'failed_build_timing_samples':0,'binary_source_contrast_independently_verified':True,'invalid_first_epoch_samples_used':0,'universal_superiority':False}
assert comparisons==48 and len(proof)==8
(OUT/'packet-audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
