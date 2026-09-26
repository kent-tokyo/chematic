"""Reject malformed annotated-tag output in release metadata timestamps."""

from scripts.generate_release_metadata import valid_released_at


def test_tag_datetime_is_accepted():
    assert valid_released_at("2026-09-26T20:27:45+09:00")
    assert valid_released_at("2026-09-26T11:27:45Z")


def test_annotated_tag_header_is_rejected():
    assert not valid_released_at("tag v1.0.27\nTagger: kent-tokyo\n2026-09-26T20:27:45+09:00")
    assert not valid_released_at("2026-09-26T99:27:45+09:00")
