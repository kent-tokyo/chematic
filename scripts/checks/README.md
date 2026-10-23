# Repository checks

This directory contains small, deterministic checks for repository contracts.
They may be run directly from the repository root and are also imported by the
script test suite.

- General API, release, naming, documentation, and workflow checks live here.
- Checks tied to a recorded measurement live in [`evidence/`](evidence/).
- Chemistry validation runners remain in [`../validation/`](../validation/).

New checks should be placed here instead of `scripts/` root. When moving an
existing check, update workflows, tests, public commands, and its repository
root lookup together.
