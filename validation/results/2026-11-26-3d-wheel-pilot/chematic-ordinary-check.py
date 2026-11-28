import json,sys,cmath
import chematic
out=[]
for row in json.load(open('/private/tmp/chematic-ordinary20-input.json')):
 try:
  m=chematic.from_smiles(row['smiles']);coords=m.generate_3d()
  out.append({'index':row['index'],'coordinates':coords})
 except Exception as e:
  out.append({'index':row['index'],'error':repr(e)})
json.dump(out,open(sys.argv[1],'w'),allow_nan=False)
