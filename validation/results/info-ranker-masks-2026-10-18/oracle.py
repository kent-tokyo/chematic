"""Seeded mask/vote schedules, live RDKit comparisons and isolated diagnostics."""
from pathlib import Path
import json
import math
from fractions import Fraction
import os
import subprocess
import sys
import numpy as np
import rdkit
from rdkit import DataStructs
from rdkit.ML.InfoTheory import rdInfoTheory as Reference
from chematic.rdkit_compat import ExplicitBitVect, rdInfoTheory as Candidate

root = Path(__file__).resolve().parent
rng = np.random.default_rng(960118)
summary = dict(rdkit_version=rdkit.__version__, seed=960118, datasets=96,
               metric_source_runs=0, matrix_comparisons=0, ranked_rows=0, independent_checks=0, empty_shape_checks=0, reference_anomalies=0,
               max_score_absolute_error=0.0, count_mismatches=0)

class Sparse:
    def __init__(self, width, bits):
        self.width, self.bits = width, tuple(bits)
    def GetNumBits(self):
        return self.width
    def GetOnBits(self):
        return iter(self.bits)

def reference_vector(width, bits, sparse):
    fp = (DataStructs.SparseBitVect if sparse else DataStructs.ExplicitBitVect)(width)
    for bit in bits:
        fp.SetBit(bit)
    return fp

def candidate_vector(width, bits, kind):
    if kind == 1:
        return Sparse(width, bits)
    if kind in (2, 3):
        return reference_vector(width, bits, kind == 3)
    fp = ExplicitBitVect(width)
    for bit in bits:
        fp.SetBit(bit)
    return fp

anomalies = []

def check_independent(candidate, eligible, model, totals, metric):
    got = {int(row[0]): row[1:] for row in candidate.GetTopN(len(eligible))}
    assert sorted(got) == sorted(eligible)
    n = sum(totals)
    def entropy(values):
        total = sum(values)
        return -math.fsum((v / total) * math.log2(v / total) for v in values if v) if total else 0.0
    for bit in eligible:
        on = model[bit]
        off = [total - hit for total, hit in zip(totals,on)]
        np.testing.assert_array_equal(got[bit][1:],on)
        if metric == 1:
            expected = entropy(totals) - sum(on) / n * entropy(on) - sum(off) / n * entropy(off)
        else:
            score = Fraction(0)
            for observed in (on, off):
                for hit, total in zip(observed, totals):
                    expected_count = Fraction(sum(observed) * total, n)
                    if expected_count:
                        score += (hit - expected_count) ** 2 / expected_count
            expected = float(score)
        assert math.isclose(float(got[bit][0]),expected,abs_tol=4e-14,rel_tol=2e-13), (got[bit][0],expected)
    summary['independent_checks'] += 1

def compare(candidate, reference, eligible, context):
    actual = candidate.GetTopN(len(eligible))
    expected = reference.GetTopN(len(eligible))
    actual = actual[np.argsort(actual[:, 0])]
    expected = expected[np.argsort(expected[:, 0])]
    np.testing.assert_array_equal(actual[:, 0], sorted(eligible))
    if expected.shape != actual.shape or list(expected[:,0]) != sorted(eligible):
        anomalies.append(dict(context=context, symptom='reference returned missing/duplicate eligible bit IDs', candidate=actual.tolist(), reference=expected.tolist()))
        summary['reference_anomalies'] += 1
        return
    np.testing.assert_array_equal(actual[:, 0], expected[:, 0])
    np.testing.assert_array_equal(actual[:, 2:], expected[:, 2:])
    np.testing.assert_allclose(actual[:, 1], expected[:, 1], atol=4e-14, rtol=2e-13)
    summary['max_score_absolute_error'] = max(summary['max_score_absolute_error'], float(np.max(np.abs(actual[:, 1] - expected[:, 1]))))
    summary['matrix_comparisons'] += 1
    summary['ranked_rows'] += len(eligible)

for dataset in range(summary['datasets']):
    width, classes = int(rng.integers(4, 33)), int(rng.integers(2, 5))
    phases = []
    for phase in range(4):
        mask = sorted(int(i) for i in rng.choice(width, size=int(rng.integers(1, width + 1)), replace=False))
        if phase == 1:
            mask = []  # Suppressed on-counts still increment class totals.
        if phase == 3:
            mask = list(range(width))
        votes = []
        for record in range(classes * 2):
            bits = [i for i in range(width) if rng.random() < 0.4]
            votes.append((bits, record % classes))
        phases.append((mask, votes))
    for metric in (1, 3):
        for kind in range(4):
            actual = Candidate.InfoBitRanker(width, classes, metric)
            expected = Reference.InfoBitRanker(width, classes, Reference.InfoType.ENTROPY if metric == 1 else Reference.InfoType.CHISQUARE)
            model = [[0] * classes for _ in range(width)]
            totals = [0] * classes
            for phase, (mask, votes) in enumerate(phases):
                # Duplicate IDs normalize to a set; bounded raw input.
                supplied = mask + mask[:1] if mask and len(mask) < width else mask
                actual.SetMaskBits(iter(supplied))
                expected.SetMaskBits(supplied)
                for bits, label in votes:
                    totals[label] += 1
                    for bit in bits:
                        if bit in mask:
                            model[bit][label] += 1
                    actual.AccumulateVotes(candidate_vector(width, bits, kind), label)
                    expected.AccumulateVotes(reference_vector(width, bits, kind % 2 == 1), label)
                if mask:
                    check_independent(actual,mask,model,totals,metric)
                    compare(actual, expected, mask, dict(dataset=dataset, width=width, classes=classes, metric=metric, source_kind=kind, phase=phase, phases=phases[:phase+1]))
                else:
                    assert actual.GetTopN(0).shape == expected.GetTopN(0).shape
                    summary['empty_shape_checks'] += 1
            summary['metric_source_runs'] += 1

# Failures are isolated from ordinary comparisons and never counted as parity.
common = '''import resource
resource.setrlimit(resource.RLIMIT_CORE, (0,0))
import json
from rdkit import DataStructs
from rdkit.ML.InfoTheory import rdInfoTheory as I
r=I.InfoBitRanker(4,2)
b=DataStructs.ExplicitBitVect(4)
b.SetBit(0)
r.AccumulateVotes(b,0)
r.SetMaskBits([0,1])
'''
diagnostics = []
for name, code in [
    ('zero_top_with_nonempty_mask', 'print(r.GetTopN(0).tolist(),flush=True)'),
    ('invalid_mask_can_replace_previous_mask', '''try:r.SetMaskBits([2,9])
except Exception as e:print(type(e).__name__,flush=True)
try:print(json.dumps(r.GetTopN(1).tolist()),flush=True)
except Exception as e:print(type(e).__name__,flush=True)'''),
]:
    run = subprocess.run([sys.executable, '-u', '-c', common + code],
                         env=dict(os.environ), capture_output=True, text=True, timeout=15)
    diagnostics.append(dict(name=name, returncode=run.returncode, stdout=run.stdout, stderr=run.stderr))
actual = Candidate.InfoBitRanker(4,2)
actual.AccumulateVotes(Sparse(4,[0]),0)
actual.SetMaskBits([0,1])
assert actual.GetTopN(0).shape == (0,4)
before = actual.GetTopN(2)
try:
    actual.SetMaskBits([2,9])
except ValueError:
    pass
else:
    raise AssertionError('invalid mask must fail')
np.testing.assert_array_equal(actual.GetTopN(2), before)
unique = {}
for anomaly in anomalies:
    context = anomaly['context']
    key = (context['dataset'], context['metric'], context['phase'])
    kind = context.pop('source_kind')
    context.pop('phases')  # Reproduced exactly from the committed seed/script.
    if key not in unique:
        unique[key] = dict(anomaly, source_kinds=[])
    unique[key]['source_kinds'].append(kind)
summary['distinct_reference_anomalies'] = len(unique)
(root / 'reference-anomalies.json').write_text(json.dumps(list(unique.values()), indent=2) + '\n')
(root / 'separate-diagnostics.json').write_text(json.dumps(diagnostics, indent=2) + '\n')
(root / 'oracle-summary.json').write_text(json.dumps(summary, indent=2) + '\n')
print(json.dumps(summary, indent=2))
print('Separate diagnostics:', [(d['name'], d['returncode']) for d in diagnostics])
