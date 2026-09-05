# Release metadata

Each release publishes a small JSON document at
`release-metadata/v<version>.json` and attaches the same file to the GitHub
Release. The schema is versioned at
[`docs/release-metadata-schema.json`](release-metadata-schema.json).

The document is intended for downstream sites and integrations that should not
scrape rendered pages. It contains release identity, registry URLs, MCP
capabilities, WASM size evidence, and links to historical benchmark records.
Benchmark entries are explicitly marked `historical`; they are not current
performance claims and must retain their pinned version, corpus, hardware, and
source path.

The release workflow regenerates the attached asset from the tag commit, so the
commit and release timestamp cannot silently drift from the tag. The checked-in
file is the stable raw GitHub path for the current public release.

Validate the checked-in document locally:

```bash
python3 scripts/check_release_metadata.py
```

Operation-level comparison scorecards are separately validated by the
dependency-free `scripts/validate_scorecard.py`; it rejects stale target
versions, missing corpus/engine provenance, and claims based on unsupported,
failed, missing, or not-measured rows.

Generate a release asset manually (the workflow supplies these values from the
tag):

```bash
python3 scripts/generate_release_metadata.py \
  --version 1.0.7 \
  --commit "$(git rev-list -1 v1.0.7)" \
  --released-at "$(git show -s --format=%cI v1.0.7)" \
  --output /tmp/chematic-release-metadata.json
```
