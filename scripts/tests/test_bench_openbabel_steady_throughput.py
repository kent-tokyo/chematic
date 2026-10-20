from pathlib import Path
import importlib.util
import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("steady_native", ROOT / "benchmarks/bench_openbabel_steady_throughput.py")
steady = importlib.util.module_from_spec(spec)
spec.loader.exec_module(steady)

def row():
    return dict(mode="write", repeats=20, records=20000, output_units=2000,
                elapsed_ns=123456789, warmup_repeats=100, warmup_elapsed_ns=510000000,
                minimum_warmup_ns=500000000, minimum_timed_ns=100000000)

def test_actual_counts_and_raw_time_survive_normalization():
    original = row()
    normalized = steady.normalize(original, "write", 3, 1000)
    assert normalized["elapsed_ns"] == 123456789 // 20
    assert normalized["elapsed_ns_raw"] == 123456789
    assert normalized["repeats"] == 20 and normalized["records"] == 20000
    assert normalized["output_units"] == 2000
    assert original == row()

@pytest.mark.parametrize("key,value", [
    ("mode", "parse"), ("repeats", 2), ("repeats", True), ("repeats", 20.0),
    ("records", 20001), ("records", True), ("elapsed_ns", 0), ("elapsed_ns", True),
    ("elapsed_ns", 99999999), ("warmup_repeats", 0), ("warmup_elapsed_ns", 499999999),
    ("minimum_warmup_ns", 0), ("minimum_timed_ns", 0),
])
def test_incomplete_or_miscounted_work_is_never_ranked(key, value):
    invalid = row(); invalid[key] = value
    with pytest.raises(ValueError): steady.normalize(invalid, "write", 3, 1000)
