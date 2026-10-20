# Fixture and report generators

`generate_*_boundary_fixtures.py` and the `gen_*` and `build_*` entry points
write committed reference fixtures, corpora, or reports. Use the exact input
and dependency versions noted in the corresponding validation or benchmark
record. Run generators from the repository root and review output diffs
before committing them. The directory moves did not regenerate any fixture.
