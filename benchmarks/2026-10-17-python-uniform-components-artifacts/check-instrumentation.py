from pathlib import Path
from types import SimpleNamespace
from collections import Counter
import importlib.util,json
O=Path(__file__).resolve().parent;R=Path('/private/tmp/chematic-python-asa-preparation-20261017')
spec=importlib.util.spec_from_file_location('fixed_batches',O/'bench-fixed-batches.py');bench=importlib.util.module_from_spec(spec);spec.loader.exec_module(bench)
next_uid=0;calls=Counter()
class Molecule:
 def __init__(self,uid):self.uid=uid
def parse(_):
 global next_uid
 next_uid+=1;return Molecule(next_uid)
def evaluate(mol):
 calls[mol.uid]+=1;return 1
bench.table=lambda engine:(parse,{'labute_asa':evaluate})
result=bench.worker(SimpleNamespace(worker='baseline',corpus=R/'validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi',offset=0,limit=1000,ops='labute_asa'))
assert next_uid==61001,next_uid
assert Counter(calls.values())=={1:40001,2:20000},Counter(calls.values())
assert len(result['operations']['labute_asa']['values'])==20000
(O/'instrumentation-check.json').write_text(json.dumps({'passed':True,'scope':'Synthetic factory/state spy; all preflight elapsed fields discarded. Not formal samples.','factory_calls':next_uid,'evaluated_once_objects':40001,'evaluated_twice_objects':20000,'returned_rows_each_mode':20000,'first_use_is_fresh':True,'repeated_use_is_warmed_once':True,'parse_inclusive_calls_factory_each_row':True},indent=2)+'\n')
print('Independent factory/state counts verified')
