from pathlib import Path
import subprocess,hashlib,json
r=Path('/private/tmp/chematic-parity-lazy-20261015');s=Path('/private/tmp/chematic-parity-lazy-source-20261015')
old=(r/'baseline-model.rs').read_text();start=old.index('pub(crate) fn build_molecule_from_model(');end=old.index('\n#[cfg(test)]',start);body=old[start:end]
expected='mod frozen {\nuse super::*;\n'+body+'\n}\n'
t=(s/'crates/chematic-perception/src/aromaticity/indexed_normalization_tests.rs').read_text();actual=t[t.index('mod frozen {'):t.index('fn compare_graph(')]
def formatted(x):return subprocess.check_output(['rustfmt','--edition','2024','--emit','stdout'],input=x,text=True)
a=formatted(expected);b=formatted(actual);assert a==b
(r/'frozen-oracle-proof.json').write_text(json.dumps({'complete_immediate_baseline_builder_order_and_finish_bodies':True,'only_rustfmt_whitespace_difference':True,'formatted_sha256':hashlib.sha256(a.encode()).hexdigest(),'baseline_body_sha256':hashlib.sha256(body.encode()).hexdigest(),'uses_new_production_table_or_finish_implementation':False},indent=2)+'\n')
print('Complete frozen baseline oracle matches after rustfmt')
