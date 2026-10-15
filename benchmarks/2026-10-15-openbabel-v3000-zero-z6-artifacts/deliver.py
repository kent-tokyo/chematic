from pathlib import Path
import json,hashlib,shutil,os
r=Path('/private/tmp/chematic-openbabel-zero-z6-20261015');src=r/'delivery';dst=Path('/Users/k_tanabe/.codex/worktrees/relevant-cycle-speed/chematic');m=json.loads((r/'source-manifest.json').read_text());index=json.loads((r/'delivery-sha256.json').read_text());sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();adopt=json.loads((r/'decision.json').read_text())['adopt_source'];changed=m['only_changed_existing_source'] if adopt else [];new=m['added_source_files']+m['added_test_files'] if adopt else []
for n,h in m['baseline_source_files_sha256'].items():
 expected={h,index[n]} if n in changed else {h}
 assert sha(dst/n) in expected,'Attached baseline changed: '+n
for n in index:
 assert n in changed+new or n=='benchmarks/2026-10-15-openbabel-v3000-zero-z6.md' or n.startswith('benchmarks/2026-10-15-openbabel-v3000-zero-z6-artifacts/'),n
 assert sha(src/n)==index[n],n
 target=dst/n
 assert not target.is_symlink(),n
 if target.exists() and n not in changed:assert sha(target)==index[n],'New destination differs: '+n
for n,h in index.items():
 target=dst/n;target.parent.mkdir(parents=True,exist_ok=True);tmp=target.with_name(target.name+'.direct-lines-copy');mode=(target if target.exists() else src/n).stat().st_mode & 0o777;shutil.copyfile(src/n,tmp);tmp.chmod(mode);assert sha(tmp)==h,n;os.replace(tmp,target);assert sha(target)==h,n
for n,h in m['baseline_source_files_sha256'].items():
 if n not in changed:assert sha(dst/n)==h,'Unrelated attached file changed: '+n
(r/'delivery-verification.json').write_text(json.dumps({'delivered_sha256':index,'unchanged_baseline_files':len(m['baseline_source_files_sha256'])-len(changed),'unchanged_previous_candidate_files_preserved':True,'no_shared_git_mutation':True,'merged_or_published':False},indent=2)+'\n')
print('Delivered',len(index),'files; all',len(m['baseline_source_files_sha256'])-len(changed),'unrelated frozen source/input files unchanged; no Git mutation')
