#!/usr/bin/env python3
"""Extract latest release notes from CHANGELOG.md and print a gh release create command."""

import re
import subprocess
import sys
from pathlib import Path

REPO = "kent-tokyo/chematic"
CHANGELOG = Path(__file__).parent.parent / "CHANGELOG.md"


def extract_latest_block(text: str) -> tuple[str, str]:
    """Return (version, body) for the first versioned section in CHANGELOG."""
    pattern = re.compile(r"^## \[(\d+\.\d+\.\d+)\]", re.MULTILINE)
    matches = list(pattern.finditer(text))
    if not matches:
        sys.exit("No versioned section found in CHANGELOG.md")
    start = matches[0]
    version = start.group(1)
    body_start = start.end()
    body_end = matches[1].start() if len(matches) > 1 else len(text)
    body = text[body_start:body_end].strip()
    return version, body


def main() -> None:
    text = CHANGELOG.read_text()
    version, body = extract_latest_block(text)
    tag = f"v{version}"

    footer = f"""
---

**Packages**
- [crates.io](https://crates.io/crates/chematic/{version}) · [docs.rs](https://docs.rs/chematic/{version})
- [PyPI](https://pypi.org/project/chematic/{version}/) · [npm](https://www.npmjs.com/package/@kent-tokyo/chematic/v/{version})
- [Live demo](https://kent-tokyo.github.io/chematic/playground/) · [Docs](https://kent-tokyo.github.io/chematic/)
""".strip()

    notes = f"{body}\n\n{footer}"

    print(f"# Release {tag}\n")
    print("Run this command to create the GitHub Release:\n")
    # Write notes to a temp file to avoid shell quoting issues
    notes_file = Path("/tmp/chematic_release_notes.md")
    notes_file.write_text(notes)
    print(f"gh release create {tag} --title '{tag}' --notes-file /tmp/chematic_release_notes.md")
    print()
    print("--- Preview ---")
    print(notes)

    if "--create" in sys.argv:
        subprocess.run(
            ["gh", "release", "create", tag, "--title", tag, "--notes-file", str(notes_file)],
            check=True,
        )
        print(f"\nRelease {tag} created.")


if __name__ == "__main__":
    main()
