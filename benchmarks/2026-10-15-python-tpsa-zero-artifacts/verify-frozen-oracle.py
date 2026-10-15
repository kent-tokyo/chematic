from pathlib import Path
import subprocess,hashlib,json
r=Path('/private/tmp/chematic-tpsa-zero-20261015');s=Path('/private/tmp/chematic-tpsa-zero-source-20261015');old=(r/'baseline-tpsa.rs').read_text();a=old[old.index('fn rdkit_tpsa_impl('):old.index('\n#[cfg(test)]',old.index('fn rdkit_tpsa_impl('))].replace('fn rdkit_tpsa_impl(', 'fn frozen_default_tpsa(');t=(s/'crates/chematic-chem/src/descriptors/tpsa_zero_tests.rs').read_text();b=t[t.index('fn frozen_default_tpsa('):t.index('\nfn compare(')]
def formatted(x):return subprocess.check_output(['rustfmt','--edition','2024','--emit','stdout'],input=x,text=True)
x=formatted(a);y=formatted(b);assert x==y;(r/'frozen-oracle-proof.json').write_text(json.dumps({'complete_immediate_baseline_default_tpsa_body':True,'only_function_name_and_rustfmt_difference':True,'formatted_sha256':hashlib.sha256(x.encode()).hexdigest(),'calls_new_guard':False},indent=2)+'\n');print('Original default TPSA body retained exactly apart from name and rustfmt')
