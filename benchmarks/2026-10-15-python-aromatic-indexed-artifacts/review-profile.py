from pathlib import Path
import re,json
r=Path('/private/tmp/chematic-parity-lazy-20261015'); result={}
terms=['descriptor_aromaticity','build_molecule_from_model','finish_molecule_with_flags','compute_bond_orders','sssr_rings_in_components','horton_candidates_fast','rdkit_model_may_disagree','assign_from_kekulized','get_atom_electron_donor_type','implicit_hcount']
for p in sorted(r.glob('sample-*.txt')):
 stack=[];den=0;counts={k:0 for k in terms}
 for line in p.read_text().splitlines():
  m=re.match(r'^([ +!|:]*)(\d+) (.*?)(?:  \(in |$)',line)
  if not m:continue
  depth=len(m[1]);count=int(m[2]);name=m[3]
  while stack and stack[-1][0]>=depth:stack.pop()
  names=[x[1] for x in stack]
  getter=lambda s:'___pymethod_get_rdkit_tpsa__' in s or '___pymethod_get_labute_asa__' in s
  if getter(name) and not any(getter(a) for a in names):den+=count
  if any(getter(a) for a in names):
   for term in terms:
    if term in name and not any(term in a for a in names):counts[term]+=count
  stack.append((depth,name))
 assert den>0,p
 result[p.name]={'getter_samples':den,'inclusive_samples_overlapping':counts}
(r/'getter-profile-review.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
