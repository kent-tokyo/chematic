import gzip
import json

import pytest

from scripts.evidence_io import (
    collect_error,
    read_gzip_jsonl,
    read_gzip_jsonl_with_sha,
    read_json,
    require,
    sha256_bytes,
)


def test_json_and_gzip_jsonl_helpers_preserve_bytes(tmp_path):
    document = tmp_path / "document.json"
    document.write_text('{"value": 3}\n', encoding="utf-8")
    assert read_json(document) == {"value": 3}

    body = b'{"input_index":0}\n{"input_index":1}\n'
    rows = tmp_path / "rows.jsonl.gz"
    rows.write_bytes(gzip.compress(body, mtime=0))
    raw, decoded = read_gzip_jsonl(rows)
    assert raw == body
    assert decoded == [{"input_index": 0}, {"input_index": 1}]
    assert read_gzip_jsonl_with_sha(rows) == (decoded, sha256_bytes(body))


def test_fail_fast_and_collecting_requirements_remain_distinct():
    require(True, "unused")
    with pytest.raises(ValueError, match="failed gate"):
        require(False, "failed gate")

    errors: list[str] = []
    collect_error(True, "unused", errors)
    collect_error(False, "collected gate", errors)
    assert errors == ["collected gate"]
