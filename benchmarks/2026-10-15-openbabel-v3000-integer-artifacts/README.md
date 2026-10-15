# V3000 integer writer evidence

Read decision.json, measurement-plan.json, source-manifest.json, all paired21 reports, continued21.json and adaptive-selection.json. Focused reports remain separate diagnostics. All original inputs and semantic exclusions are retained.

Extract candidate-build-source.tar.gz: candidate/ contains the full build/input source, baseline-overrides/ contains original changed files, corpora/ contains fixed exact benchmark streams. Reconstruct baseline by applying overrides and deleting added files listed in the manifest. verification.json checks all source hashes. candidate.patch reconstructs the candidate in an independent tiny Git repository. No new molecular corpus rows are added; original corpus provenance remains in chemistry-speed-followups-v1.manifest.json.

The two native macOS arm64 executables are exact measured binaries. Rebuild with the recorded Rust toolchain and cargo build --release --locked -p chematic-mol --example openbabel_v3000_throughput, release opt3/LTO/codegen-units1. Bit identity of a rebuild is not promised. Open Babel native harness/library/plugin paths and hashes are pinned; its C++ harness source is archived.

Use extracted scripts/bench_openbabel_v3000_throughput.py with --baseline, --candidate, --openbabel, --corpus, --records 1000, --blocks 21, --repeats 5, --source-manifest and --output in the pinned observer environment. Use --records 64 for controlled chains. Darwin time -l requires OS-statistics access. Driver paths are host-specific. Gzip streams are verified against every raw file. artifact-sha256.json indexes other packet files. These are unmerged version 1.1.4 source candidates.
