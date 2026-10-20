"""Published wheels must cover the declared Python and platform matrix."""

import unittest

from scripts.release.check_pypi_dist import PLATFORMS, PYTHONS, check_layout


class DistributionTests(unittest.TestCase):
    def setUp(self):
        self.names = {"chematic-1.2.3.tar.gz"}
        platform_tags = {
            "linux-x86_64": "manylinux_2_17_x86_64.manylinux2014_x86_64",
            "linux-aarch64": "manylinux_2_17_aarch64.manylinux2014_aarch64",
            "macos-x86_64": "macosx_10_12_x86_64",
            "macos-arm64": "macosx_11_0_arm64",
            "windows-amd64": "win_amd64",
        }
        for python in PYTHONS:
            for platform in PLATFORMS:
                self.names.add(f"chematic-1.2.3-cp{python}-cp{python}-{platform_tags[platform]}.whl")

    def test_complete_matrix(self):
        check_layout(self.names, "1.2.3")

    def test_missing_wheel_fails(self):
        self.names.remove("chematic-1.2.3-cp39-cp39-win_amd64.whl")
        with self.assertRaisesRegex(ValueError, "missing"):
            check_layout(self.names, "1.2.3")

    def test_wrong_release_fails(self):
        with self.assertRaises(ValueError):
            check_layout(self.names, "1.2.4")


if __name__ == "__main__":
    unittest.main()
