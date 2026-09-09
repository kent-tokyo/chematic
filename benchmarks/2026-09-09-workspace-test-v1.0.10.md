# Workspace test gate: v1.0.10 candidate

After reclaiming only regenerable Cargo build output, the complete Rust
workspace was tested offline with bounded build parallelism.

```sh
CARGO_BUILD_JOBS=2 cargo test --workspace --offline
```

The run completed successfully across workspace unit tests, integration tests,
and doctests. Aggregating the 94 `test result` lines produced:

| Result | Count |
|---|---:|
| passed | 4,849 |
| failed | 0 |
| ignored | 34 |
| measured | 0 |
| filtered out | 0 |

The ignored tests are deliberate long-running or environment-specific lanes;
they are not converted into passes. This is local Rust evidence and does not
cover browser artifact regeneration, external review, registry publication, or
future release actions.

Machine-readable record: [`2026-09-09-workspace-test-v1.0.10.json`](2026-09-09-workspace-test-v1.0.10.json).
