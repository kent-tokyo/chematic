import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

SPEC = importlib.util.spec_from_file_location(
    "morgan_linear", Path(__file__).parents[1] / "benchmarks" / "bench_python_morgan_linear.py"
)
bench = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(bench)


class MorganBatchTests(unittest.TestCase):
    def measured(self, mode):
        objects, calls, events = [], [], []

        def parse(s):
            obj = {"row": s}
            objects.append(obj)
            events.append("parse")
            return obj

        def fn(obj):
            calls.append(obj)
            events.append("fn")
            return obj["row"]

        def clock():
            events.append("clock")
            return len([x for x in events if x == "clock"]) * 10

        record, values = bench.measure_mode(
            parse, fn, lambda x: x, ["C", "O"], mode, 25, clock=clock
        )
        self.assertEqual(values, ["C", "O"])
        self.assertEqual(record["batch_ns"], [10, 10, 10])
        self.assertEqual(record["processed_rows"], 6)
        self.assertEqual(record["elapsed_ns_per_batch"], 10)
        return objects, calls, events

    def test_first_use_parses_fresh_outside_every_timer(self):
        objects, calls, events = self.measured("first_use")
        self.assertEqual(len(objects), 6)
        self.assertTrue(all(a is b for a, b in zip(objects, calls)))
        self.assertEqual(events, ["parse", "parse", "clock", "fn", "fn", "clock"] * 3)

    def test_repeated_reuses_objects_after_untimed_priming(self):
        objects, calls, events = self.measured("repeated_use")
        self.assertEqual(len(objects), 2)
        self.assertTrue(all(obj is objects[i % 2] for i, obj in enumerate(calls)))
        self.assertEqual(
            events, ["parse", "parse", "fn", "fn"] + ["clock", "fn", "fn", "clock"] * 3
        )

    def test_inclusive_parses_inside_every_timer(self):
        objects, calls, events = self.measured("parse_inclusive")
        self.assertEqual(len(objects), 6)
        self.assertEqual(events, ["clock", "parse", "fn", "parse", "fn", "clock"] * 3)

    def test_normalizes_actual_batches_with_fractional_remainder(self):
        ticks = iter([0, 7, 10, 20])
        record, _ = bench.measure_mode(
            str, str, str, ["C"], "first_use", 11, clock=lambda: next(ticks)
        )
        self.assertEqual(record["elapsed_ns_raw"], 17)
        self.assertEqual(record["elapsed_ns_per_batch"], 8)
        self.assertEqual(record["batches"], 2)

    def test_refusals_nondeterminism_and_nonpositive_time_fail(self):
        for parse, fn, clock in [
            (lambda _: None, str, iter([0, 1])),
            (str, lambda _: "changed", iter([0, 1])),
            (str, str, iter([0, 0])),
        ]:
            with self.assertRaises(ValueError):
                bench.measure_mode(
                    parse,
                    fn,
                    str,
                    ["C"],
                    "first_use",
                    1,
                    expected=["C"],
                    clock=lambda: next(clock),
                )

    def test_accounting_corruption_is_rejected(self):
        ticks = iter([0, 10])
        record, _ = bench.measure_mode(
            str, str, str, ["C"], "first_use", 1, clock=lambda: next(ticks)
        )
        for field in [
            "elapsed_ns_raw",
            "processed_rows",
            "batches",
            "elapsed_ns_per_batch",
            "minimum_accumulated_ns",
        ]:
            bad = dict(record)
            bad[field] += 1
            with self.assertRaises(ValueError):
                bench.check_timing(bad, "first_use", 1, 1)
        bad = dict(record, output_sha256_per_batch=[])
        with self.assertRaises(ValueError):
            bench.check_timing(bad, "first_use", 1, 1)

    def test_unfiltered_window_and_atomic_failure_preserve_evidence(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "rows"
            path.write_text("C first\ninvalid second\nO third\n")
            self.assertEqual(bench.rows_of(path, 1, 2), ["invalid", "O"])
            with self.assertRaises(ValueError):
                bench.rows_of(path, 2, 2)
            result = Path(directory) / "journal.json"
            bench.atomic_json(result, {"original": True})
            old = result.read_bytes()
            with patch.object(
                Path, "replace", side_effect=OSError("simulated failure")
            ):
                with self.assertRaises(OSError):
                    bench.atomic_json(result, {"new": True})
            self.assertEqual(result.read_bytes(), old)


if __name__ == "__main__":
    unittest.main()
