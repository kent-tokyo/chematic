"""The CI shard planner must include new integration tests automatically."""

import unittest

from scripts.checks.ci_integration_shard import assign_packages


class ShardTests(unittest.TestCase):
    def test_all_test_packages_are_assigned_once(self):
        metadata = {
            "workspace_members": ["core", "mol", "new", "empty"],
            "packages": [
                {"id": "core", "name": "core", "targets": [{"kind": ["test"]}]},
                {"id": "mol", "name": "mol", "targets": [{"kind": ["test"]}] * 3},
                {"id": "new", "name": "new", "targets": [{"kind": ["test"]}]},
                {"id": "empty", "name": "empty", "targets": [{"kind": ["lib"]}]},
                {"id": "external", "name": "external", "targets": [{"kind": ["test"]}]},
            ],
        }
        shards = assign_packages(metadata, 2)
        self.assertEqual(sorted(sum(shards, [])), ["core", "mol", "new"])
        self.assertEqual(shards, assign_packages(metadata, 2))
        self.assertTrue(all(shards))

    def test_rejects_empty_shards(self):
        metadata = {"workspace_members": [], "packages": []}
        with self.assertRaises(ValueError):
            assign_packages(metadata, 2)


if __name__ == "__main__":
    unittest.main()
