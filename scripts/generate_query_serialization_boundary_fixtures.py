"""Freeze native query/reaction SMARTS serialization and parse rejection."""
import itertools,json
from pathlib import Path
from rdkit import Chem,rdBase
from rdkit.Chem import rdChemReactions
assert rdBase.rdkitVersion=='2026.03.1'
queries=['B','P','S','I','b1ccccc1','o1cccc1','s1cccc1','p1ccccc1','[a]','[A]',
'[13]','[2H]','[13H]','[H:1]','[C++]','[O--]','[N-3]','[C+3]','[d2]','[r]','[k]',
'[D{2-4}]','[D{2-}]','[D{-4}]','[!D{2-4}]','[r5,r6]','[!C;!N]','[!$(C=O)]',
'[$(C=O),$(N#N)]','[C,N;H1,H2]','[C;!$(C=O)]','[!C&!N]','[C;N,O]',
'C-&!@C','C-;!@C','C!@C','C!,=C','C!-,!#C','C~C','C/C=C\\C',
'C\\C=C/C','C/1CCCCC1','C=1CCCCC1','C%12CCCCC%12','[se]1cccc1',
'[as]1ccccc1','[te]1cccc1','[C:0]','[C:12]','C=','C%','C%1','C%a1',
'[$(CC]','[D{4-2}]','[D{a-2}]','[C','C(C','C1CC','[foo]','[C:]',
'[C;$(C(=O)O)]','[!$([N;H0]);H1]']
# Real logical combinations cross atom and bond negation/precedence, not
# manually constructed internal query trees.
for a,b,bond in itertools.product(['C','N','[a]','[A]','[13C]','[!C]'],['O','S','[N+]'],['-','=','#','-&!@','-,=','!#']):
    queries.append(a+bond+b)
rows=[]
for text in dict.fromkeys(queries):
    m=Chem.MolFromSmarts(text)
    rows.append(dict(input=text,expected=Chem.MolToSmarts(m) if m is not None else None))
reactions=['[C:1]>>[C:1]','[C:1]>O>[C:1]','[N:1].[C:2]>>[N:1][C:2]',
           '([C:1].[O:2])>>[C:1][O:2]','[C:1]/[C:2]=[C:3]\\[C:4]>>[C:1]/[C:2]=[C:3]/[C:4]',
           '>>','C>>','>>O','C>O>','[13C:1]>>[13C:1]','[a:1]>>[a:1]']
rr=[]
for text in reactions:
    try:r=rdChemReactions.ReactionFromSmarts(text);expected=rdChemReactions.ReactionToSmarts(r) if r is not None else None
    except ValueError:expected=None
    rr.append(dict(input=text,expected=expected))
(Path(__file__).resolve().parents[1]/'validation/rdkit-2026.03.1-query-serialization-boundary.json').write_text(json.dumps(dict(rdkit_version=rdBase.rdkitVersion,queries=rows,reactions=rr),indent=2)+'\n')
print(len(rows),'queries',len(rr),'reactions')
