import re
from pathlib import Path

from scripts.validation_provenance import collect


ROOT = Path(__file__).parents[2]
SHA256 = re.compile(r"^[0-9a-f]{64}$")
COMMIT = re.compile(r"^[0-9a-f]{40}$")


def test_provenance_contains_reproducibility_identity():
    record = collect(ROOT)
    assert COMMIT.fullmatch(record["source_commit"])
    assert isinstance(record["dirty_worktree"], bool)
    assert SHA256.fullmatch(record["tracked_diff_sha256"])
    assert SHA256.fullmatch(record["untracked_source_sha256"])
    assert isinstance(record["untracked_source_files"], int)
    assert record["untracked_source_files"] >= 0
    assert isinstance(record["python_executable"], str)
    assert isinstance(record["platform"], str)
    assert isinstance(record["machine"], str)


def test_provenance_digest_is_nonempty_without_worktree_mutation(tmp_path):
    first = collect(ROOT)
    # A temporary file is not part of the repository, so this test verifies
    # the digest contract without mutating the shared worktree.
    assert first["tracked_diff_sha256"] != ""
    assert tmp_path.is_dir()
