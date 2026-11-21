"""Partition integrity checks for audited COD holdouts."""
import importlib.util
import hashlib
import json
from pathlib import Path
import unittest


MODULE = Path(__file__).resolve().parents[1] / "validation" / "cod_holdout_shards.py"
spec = importlib.util.spec_from_file_location("cod_holdout_shards", MODULE)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
merge_spec = importlib.util.spec_from_file_location(
    "merge_cod_holdout_shards", MODULE.with_name("merge_cod_holdout_shards.py")
)
merge_mod = importlib.util.module_from_spec(merge_spec)
merge_spec.loader.exec_module(merge_mod)


class ShardTest(unittest.TestCase):
    def setUp(self):
        self.split = {
            "source_sha256": "source",
            "invalid_exclusions": [],
            "eligibility": {"eligible_holdout_records": 6},
            "rows": [
                {"index": i, "name": f"record-{i}", "split": "holdout", "scaffold_sha256": f"{i // 2:064x}"}
                for i in range(6)
            ],
        }
        self.split_bytes = mod.json_bytes(self.split)
        self.audit = {
            "split_sha256": mod.digest(self.split_bytes),
            "source_sha256": "source",
            "invalid_exclusion_count": 0,
            "overlaps": [],
            "rule": "parent clean",
        }
        self.audit_bytes = mod.json_bytes(self.audit)

    def test_partitions_are_disjoint_complete_and_scaffold_intact(self):
        outputs = mod.shard(self.split_bytes, self.audit_bytes, 3)
        seen = set()
        scaffold_parts = {}
        for part, (split_bytes, audit_bytes) in enumerate(outputs):
            split, audit = json.loads(split_bytes), json.loads(audit_bytes)
            self.assertEqual(audit["split_sha256"], mod.digest(split_bytes))
            rows = [r for r in split["rows"] if r["split"] == "holdout"]
            self.assertEqual(split["eligibility"]["eligible_holdout_records"], len(rows))
            for row in rows:
                self.assertNotIn(row["index"], seen)
                seen.add(row["index"])
                self.assertEqual(scaffold_parts.setdefault(row["scaffold_sha256"], part), part)
        self.assertEqual(seen, set(range(6)))

    def test_rejects_mismatched_or_unclean_audit(self):
        self.audit["overlaps"] = ["overlap"]
        with self.assertRaises(ValueError):
            mod.shard(self.split_bytes, mod.json_bytes(self.audit), 3)
        self.audit["overlaps"] = []
        self.audit["split_sha256"] = "wrong"
        with self.assertRaises(ValueError):
            mod.shard(self.split_bytes, mod.json_bytes(self.audit), 3)

    def test_merge_preserves_rows_and_rejects_missing_or_duplicate_data(self):
        shards = mod.shard(self.split_bytes, self.audit_bytes, 3)
        pairs = []
        for i, (split_bytes, _) in enumerate(shards):
            split = json.loads(split_bytes)
            rows = [dict(index=r["index"], name=r["name"], value=r["index"] * 2)
                    for r in split["rows"] if r["split"] == "holdout"]
            result = {"source_sha256": "source", "split_sha256": mod.digest(split_bytes),
                      "rdkit": "fixed", "platform": f"runner-{i}", "results": rows}
            pairs.append((split_bytes, mod.json_bytes(result)))
        merged = merge_mod.merge(self.split_bytes, pairs, "fixed")
        self.assertEqual(merged["split_sha256"], mod.digest(self.split_bytes))
        self.assertEqual([r["value"] for r in merged["results"]], [0, 2, 4, 6, 8, 10])
        self.assertEqual(len(merged["shard_provenance"]), 3)
        broken = list(pairs)
        record = json.loads(broken[0][1])
        record["results"].pop()
        broken[0] = (broken[0][0], mod.json_bytes(record))
        with self.assertRaises(ValueError):
            merge_mod.merge(self.split_bytes, broken, "fixed")
        with self.assertRaises(ValueError):
            merge_mod.merge(self.split_bytes, [pairs[0], pairs[0], pairs[2]], "fixed")

    def test_merge_rejects_different_settings_and_wrong_split_hash(self):
        pairs = []
        for split_bytes, _ in mod.shard(self.split_bytes, self.audit_bytes, 3):
            split = json.loads(split_bytes)
            rows = [dict(index=r["index"], name=r["name"]) for r in split["rows"] if r["split"] == "holdout"]
            pairs.append((split_bytes, mod.json_bytes({"source_sha256": "source", "split_sha256": mod.digest(split_bytes), "seed": 42, "platform": "linux", "results": rows})))
        changed = json.loads(pairs[1][1])
        changed["seed"] = 43
        with self.assertRaises(ValueError):
            merge_mod.merge(self.split_bytes, [pairs[0], (pairs[1][0], mod.json_bytes(changed)), pairs[2]], "fixed")
        changed["seed"] = 42
        changed["split_sha256"] = "wrong"
        with self.assertRaises(ValueError):
            merge_mod.merge(self.split_bytes, [pairs[0], (pairs[1][0], mod.json_bytes(changed)), pairs[2]], "fixed")

    def test_timed_rows_merge_with_per_shard_provenance(self):
        pairs = []
        for i, (split_bytes, _) in enumerate(mod.shard(self.split_bytes, self.audit_bytes, 3)):
            split = json.loads(split_bytes)
            rows = [dict(index=r["index"], name=r["name"]) for r in split["rows"] if r["split"] == "holdout"]
            data = {"source_sha256": "source", "seconds": 0.2, "platform": f"linux-{i}", "results": rows}
            pairs.append((split_bytes, mod.json_bytes(data)))
        merged = merge_mod.merge(self.split_bytes, pairs, "time")
        self.assertEqual(len(merged["results"]), 6)
        self.assertEqual([r["platform"] for r in merged["shard_provenance"]], ["linux-0", "linux-1", "linux-2"])
        expected = sorted(range(6), key=lambda i: (hashlib.sha256(f"record-{i}".encode()).digest(), i))
        self.assertEqual([r["index"] for r in merged["results"]], expected)


if __name__ == "__main__":
    unittest.main()
