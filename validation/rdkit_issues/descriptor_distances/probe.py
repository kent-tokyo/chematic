"""Live ordinary RDKit comparisons; independent math.dist numeric controls."""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import random

import chematic
import chematic.chematic as extension
import numpy as np
from rdkit import rdBase
from rdkit.DataManip.Metric import rdMetricMatrixCalc as metric


def expected(rows):
    return [math.dist(rows[i], rows[j]) for i in range(1,len(rows)) for j in range(i)]


def check(actual, reference):
    assert len(actual) == len(reference)
    for a,b in zip(actual, reference):
        assert math.isfinite(a) and math.isclose(a,b,rel_tol=2e-12,abs_tol=1e-12), (a,b)


def finite_json(values):
    return [float(v) if math.isfinite(v) else repr(float(v)) for v in values]


def main(output):
    records = []
    runs = 0
    for count, dimensions, scale, seed in itertools.product(
        [1,2,3,8,16], [0,1,3,16,64], [1e-3,1.0,1e3], range(4)
    ):
        rng = random.Random(seed)
        rows = [[rng.uniform(-scale,scale) for _ in range(dimensions)] for _ in range(count)]
        reference = expected(rows)
        oracle = list(metric.GetEuclideanDistMat(rows))
        check(oracle, reference)
        backing = np.zeros((count,dimensions*2),dtype=np.float64)
        backing[:,::2] = rows
        view = backing[:,::2]
        variants = [rows, tuple(map(tuple,rows)), (iter(row) for row in rows), view]
        native_results = []
        for variant in variants:
            result = chematic.euclidean_distance_matrix_checked(variant)
            check(result, reference)
            native_results.append(result)
            runs += 1
        permuted = rows[::-1]
        check(chematic.euclidean_distance_matrix_checked(permuted), expected(permuted))
        runs += 1
        records.append(dict(rows=rows, expected=reference, rdkit=oracle,
                            native=native_results[0], count=count, dimensions=dimensions,
                            scale=scale, seed=seed))

    # Exercise the intended descriptor workflow without claiming descriptor
    # agreement: both distance engines receive the same native descriptors.
    smiles = ['C','CO','CCO','CC(=O)O','c1ccccc1','CC(=O)Oc1ccccc1C(=O)O']
    molecular_rows = [[m.mw,m.logp,m.tpsa] for m in map(chematic.from_smiles,smiles)]
    native = chematic.euclidean_distance_matrix_checked(molecular_rows)
    check(native, expected(molecular_rows))
    check(list(metric.GetEuclideanDistMat(molecular_rows)), expected(molecular_rows))

    extremes = []
    for rows in [[[0],[1e200]], [[0],[1e-200]], [[0],[math.ulp(0.0)]],
                 [[0,0],[1e308,1e308]], [[0],[float.fromhex('0x1.fffffffffffffp+1023')]]]:
        reference = expected(rows)
        native = chematic.euclidean_distance_matrix_checked(rows)
        # Absolute tolerance is deliberately zero: losing a tiny norm to zero
        # must fail this check, even though its absolute error is very small.
        for a,b in zip(native,reference):
            assert a > 0 and math.isclose(a,b,rel_tol=2e-15,abs_tol=0)
        oracle = list(metric.GetEuclideanDistMat(rows))
        extremes.append(dict(rows=rows, expected=reference, native=native, rdkit=finite_json(oracle)))

    result = dict(
        rdkit_version=rdBase.rdkitVersion, numpy_version=np.__version__,
        native_path=extension.__file__,
        native_sha256=hashlib.sha256(Path(extension.__file__).read_bytes()).hexdigest(),
        ordinary_cases=len(records), ordinary_native_runs=runs,
        ordinary_native_max_absolute_error=max((abs(a-b) for r in records for a,b in zip(r['native'],r['expected'])),default=0),
        ordinary_rdkit_max_absolute_error=max((abs(a-b) for r in records for a,b in zip(r['rdkit'],r['expected'])),default=0),
        relative_tolerance=2e-12, absolute_tolerance=1e-12,
        molecular_descriptor_smiles=smiles,molecular_descriptor_rows=molecular_rows,
        molecular_descriptor_workflow_cases=1,
        extreme_cases=extremes, extreme_relative_tolerance=2e-15, extreme_absolute_tolerance=0,
        failures=0, records=records,
    )
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({key:value for key,value in result.items() if key!='records'},indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args().output)
