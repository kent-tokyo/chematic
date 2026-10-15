from pathlib import Path
import json,hashlib,shutil,importlib.util
r=Path('/private/tmp/chematic-openbabel-ez-gate-20261015');s=Path('/private/tmp/chematic-openbabel-ez-gate-source-20261015');d=r/'delivery';a=d/'benchmarks/2026-10-15-openbabel-ez-inactive-artifacts';sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();m=json.loads((r/'source-manifest.json').read_text());adopt=json.loads((r/'decision.json').read_text())['adopt_source']
for name in ['package.py','deliver.py','package.log','finalize-packet.py']:shutil.copyfile(r/name,a/name)
assert all(sha(s/n)==h for n,h in m['frozen_build_source_sha256'].items())
if adopt:
 for n in m['only_changed_existing_source']+m['added_source_files']+m['added_test_files']:assert sha(d/n)==m['frozen_build_source_sha256'][n]
index={p.relative_to(a).as_posix():sha(p) for p in sorted(a.rglob('*')) if p.is_file() and p.name!='artifact-sha256.json'};(a/'artifact-sha256.json').write_text(json.dumps(index,indent=2)+'\n');assert all(sha(a/n)==h for n,h in index.items())
outer={p.relative_to(d).as_posix():sha(p) for p in sorted(d.rglob('*')) if p.is_file()};(r/'delivery-sha256.json').write_text(json.dumps(outer,indent=2)+'\n');assert all(sha(d/n)==h for n,h in outer.items())
spec=importlib.util.spec_from_file_location('names',s/'scripts/check_product_names.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);assert not mod.scan_paths(d,[d/'benchmarks/2026-10-15-openbabel-ez-inactive.md',a/'README.md']+([d/'CHANGELOG.md'] if adopt else []))
print('Final packet verified:',len(index),'artifacts;',len(outer),'delivery files;',len(m['frozen_build_source_sha256']),'source/input hashes; adopted=',adopt)
