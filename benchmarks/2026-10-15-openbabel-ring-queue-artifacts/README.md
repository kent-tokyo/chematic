# Inline ring-search evidence

Read decision.json, measurement-plan.json, source-manifest.json, all paired21 matrices, continued21.json and adaptive-selection.json. Focused reports are separate diagnostics. Original inputs and semantic exclusions remain.

Extract candidate-build-source.tar.gz: candidate/ is the complete build/input source, baseline-overrides/ contains pre-change files, corpora/ contains exact fixed streams. Reconstruct baseline by applying overrides and removing the added test listed in the manifest. verification.json checks every source/input hash. candidate.patch reconstructs changes independently. Existing input provenance is in chemistry-speed-followups-v1.manifest.json; Main fixed corpora are unchanged; the separate synthetic workload has its own provenance manifest and plan.

The two native macOS arm64 executables are exactly measured. Rebuild with the recorded Rust toolchain and cargo build --release --locked -p chematic-mol --example openbabel_v3000_throughput, opt3/fat-LTO/codegen-units1. Bit identity of a rebuild is not promised. Open Babel harness/library/plugin hashes are pinned and C++ harness source is archived.

Run scripts/bench_openbabel_v3000_throughput.py with --baseline, --candidate, --openbabel, --corpus, --records 1000, --blocks 21, --repeats 5, --source-manifest and --output in the pinned observer environment. Controlled chains use --records 64. Darwin time -l requires OS-statistics access. Driver paths are host-specific. Exact gzip streams are verified; artifact-sha256.json indexes every other packet file. Read the adoption decision before using an archived prototype.
