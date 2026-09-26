"""Release notes must contain the entry body, not its heading/date."""

from scripts.make_release_notes import extract_latest_block


def test_extract_latest_dated_entry_without_heading_suffix():
    changelog = (
        "## [Unreleased]\n\n"
        "## [1.0.27] - 2026-09-26\n\n"
        "- Added atom tags.\n\n"
        "## [1.0.26] - 2026-09-25\n\n"
        "- Previous release.\n"
    )
    assert extract_latest_block(changelog) == ("1.0.27", "- Added atom tags.")
