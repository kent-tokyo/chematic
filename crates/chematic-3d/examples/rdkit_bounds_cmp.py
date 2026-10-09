import sys, json
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import rdDistGeom
RDLogger.DisableLog("rdApp.*")
cnt = {}
first = {}
for line in open(sys.argv[1]):
    d = json.loads(line)
    m = Chem.MolFromSmiles(d["smiles"])
    if m is None:
        continue
    mh = Chem.AddHs(m)
    for key, mac in (("std", False), ("mac", True)):
        for stage, smooth in (("raw", False), ("smooth", True)):
            k = f"{key}_{stage}"
            c = cnt.setdefault(k, [0, 0])
            c[1] += 1
            try:
                ref = rdDistGeom.GetMoleculeBoundsMatrix(mh, set15bounds=True, scaleVDW=False,
                                                         doTriangleSmoothing=smooth, useMacrocycle14config=mac)
            except Exception as e:
                ref = None
            if k not in d:
                first.setdefault(k, (d["smiles"], "missing", d.get(f"{key}_error")))
                continue
            got = np.array([float(x) for x in d[k]]).reshape(ref.shape) if ref is not None else None
            if ref is not None and np.array_equal(ref, got):
                c[0] += 1
            else:
                if ref is not None:
                    diff = np.argwhere(ref != got)
                    i, j = diff[0]
                    first.setdefault(k, (d["smiles"], int(i), int(j), ref[i, j], got[i, j], len(diff)))
                else:
                    first.setdefault(k, (d["smiles"], "rdkit failed"))
for k, v in cnt.items():
    print(k, f"{v[0]}/{v[1]}")
for k, v in first.items():
    print(k, v)
