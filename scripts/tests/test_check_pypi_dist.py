"""Published wheels must cover the declared Python and platform matrix."""

import io
import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from scripts.release.check_pypi_dist import (
    PLATFORMS,
    PYTHONS,
    check_layout,
    local_hashes,
    published_hashes,
)


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

    def test_publish_attestations_are_not_distributions(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            distribution = root / "chematic-1.2.3.tar.gz"
            distribution.write_bytes(b"sdist")
            (root / f"{distribution.name}.publish.attestation").write_text("signed")
            self.assertEqual(set(local_hashes(root)), {distribution.name})

    def test_orphaned_publish_attestation_fails(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "missing.whl.publish.attestation").write_text("signed")
            with self.assertRaisesRegex(ValueError, "orphaned"):
                local_hashes(root)

    def test_unknown_release_file_fails(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "notes.txt").write_text("unexpected")
            with self.assertRaisesRegex(ValueError, "unexpected"):
                local_hashes(root)

    def test_published_hashes_bypasses_stale_cdn_entries(self):
        payload = {
            "urls": [{"filename": "chematic.whl", "digests": {"sha256": "abc"}}]
        }
        response = io.BytesIO(json.dumps(payload).encode())
        with patch(
            "scripts.release.check_pypi_dist.urllib.request.urlopen", return_value=response
        ) as urlopen:
            self.assertEqual(published_hashes("1.2.3"), {"chematic.whl": "abc"})
        request = urlopen.call_args.args[0]
        self.assertIn("?release-verification=", request.full_url)
        self.assertEqual(request.get_header("Cache-control"), "no-cache")


if __name__ == "__main__":
    unittest.main()
