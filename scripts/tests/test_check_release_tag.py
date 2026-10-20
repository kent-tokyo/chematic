"""Release tags must identify the version checked into the same commit."""

import unittest

from scripts.check_release_tag import version_from_tag


class ReleaseTagTests(unittest.TestCase):
    def test_stable_and_prerelease(self):
        self.assertEqual(version_from_tag("v1.2.3", "1.2.3"), "1.2.3")
        self.assertEqual(version_from_tag("v1.2.3-rc.2", "1.2.3"), "1.2.3")

    def test_rejects_branch_and_wrong_version(self):
        for tag in ("main", "v1.2.3oops", "v1.2.4", "v1.2", "v1.2.3/evil"):
            with self.subTest(tag=tag), self.assertRaises(ValueError):
                version_from_tag(tag, "1.2.3")


if __name__ == "__main__":
    unittest.main()
