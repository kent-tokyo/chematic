from pathlib import Path
import gzip,hashlib,json,re,subprocess,sys,tempfile
ROOT=Path('/private/tmp/chematic-sdf-buffer-scan-20261016')
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
raw=subprocess.check_output(['git','cat-file','--batch'],input=''.join(manifest['base_revision']+':'+p+'\n' for p in paths).encode(),cwd=ROOT)
for path in paths:
 end=raw.index(b'\n',at);header=raw[at:end].split();at=end+1
 if header[-1]==b'missing':assert path in manifest['added'],path;continue
 size=int(header[2]);blob=raw[at:at+size];at+=size+1
 if path not in manifest['changed']:assert hashlib.sha256(blob).hexdigest()==manifest['source_files_sha256'][path],path
assert at==len(raw)
old_writer=subprocess.check_output(['git','show',manifest['base_revision']+':crates/chematic-mol/src/mol3000/writer.rs'],cwd=ROOT,text=True)
original=old_writer[old_writer.index('fn append_atoms_3d('):old_writer.index('fn v3000_bond_order(')].strip()
frozen=(ROOT/'crates/chematic-mol/src/mol3000/coordinate_format_tests.rs').read_text()
frozen=frozen[frozen.index('fn original_atoms('):frozen.index('fn compare(')].strip().replace('fn original_atoms(','fn append_atoms_3d(')
assert re.sub(r'\s+','',original)==re.sub(r'\s+','',frozen)
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
result={'source_file_hashes_verified':len(manifest['source_files_sha256']),'native_library_entries_verified':len(manifest['native_libraries_sha256']),'patch_reconstructed_from':manifest['base_revision'],'all_comparisons_recomputed':comparisons,'profiles':proof,'fresh_post_timing_native_observers_not_timing_samples':True,'source_post_timing_bytes_and_snapshots_equal':True,'openbabel_post_timing_declared_semantics_equal':True,'all_unnamed_source_inputs_match_base_blobs':True,'frozen_original_atom_writer_reverified_from_git':True,'common_probe_reverified_after_rustfmt':True,'baseline_has_no_ignored_python_cache_files':True,'protocol_seal_verified':True,'initial_failed_launch_completed_blocks':0,'universal_superiority':False}
assert comparisons==48 and len(proof)==8
(OUT/'packet-audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
