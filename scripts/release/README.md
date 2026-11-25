# Release checklist

The [Release workflow](../../.github/workflows/publish-pypi.yml) qualifies the
tagged source, builds and publishes packages, and creates the GitHub Release.

## Before tagging

1. Compare public Rust, Python, WASM declarations and shared JSON contracts
   with the previous published tag. For v1.2.5, use `v1.2.4`. Check
   Cleanup and the new checked conformer bindings against their documented
   contracts. See [Cleanup migration](../../docs/reference/cleanup-migration.md)
   and [checked ensembles](../../docs/reference/checked-conformer-ensembles.md).
2. Synchronize workspace and package versions, dependency pins, lockfiles,
   citation metadata and release examples. Move the relevant Unreleased notes
   into the versioned CHANGELOG section. Run the
   [release documentation check](../checks/check_release_docs_consistency.py)
   and `mkdocs build --strict`.
3. Review the versioned published-wheel chemistry expectation manifest
   required by preflight. Check the complete distribution matrix and Python 3.9
   runner support. Verify the required checks on the versioned `main` commit.
4. Tag and push the version-matching commit on `main`. An RC tag qualifies
   source without publishing registry packages.

## After tagging

Confirm source qualification, coverage and security on the tag. Verify
crates.io checksums and source commit, npm version and integrity, and the
25 PyPI wheels plus sdist against the build hashes. Check the published wheel
smoke and chemistry jobs, the GitHub Release and its metadata asset, and the
version shown on docs.rs and the project site.

If only some channels publish, check each registry before retrying. Retry
against the same tag and, where possible, the original artifacts. If the
PyPI JSON API shows the files before pip's index does, wait for index
propagation before rerunning install tests. The
[v1.2.2 channel record](../../validation/results/release-channel-verification-v1.2.2.json)
shows the evidence to retain. The v1.2.3 and v1.2.4 records still contain
pre-publication states; see the [validation report](../../docs/reference/validation.md).

Keep generated metadata, provenance, and SBOM files with the release
artifacts. Follow the [release key custody guide](../../docs/reference/release-key-custody.md)
for signing.
