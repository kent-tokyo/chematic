# Inline ring scratch evidence

Read `decision.json`, `measurement-plan.json`, `source-manifest.json`, both `frontier-*.json` reports and `adaptive-selection.json`. All blocks, modes, output gates, load observations and installed extension hashes remain. Any `focused-*.json` reports are separate adaptive diagnostics.

Extract `candidate-build-source.tar.gz` to obtain `source/`. To reconstruct the immediate baseline, replace paths using `baseline-overrides/` and remove the added test listed in the source manifest. `source-reconstruction.json` checks every source/input hash. The patch reproduces the candidate changes. Source version is 1.1.4; this is not main or a release.

The two wheel directories contain the exact measured binaries. Install each separately using Python 3.13 and `pip install --no-deps --target <directory> <wheel>`. Pin comparator dependencies from `comparator-environment.txt`. Run the extracted `scripts/bench_python_speed_frontier.py` with `--corpus source/validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi --offset 5000 --limit 1000 --blocks 20` and supply `--baseline-wheel`, `--baseline-path`, `--candidate-wheel`, `--candidate-path`, `--source-manifest`, `--base-revision` and `--output`. Repeat at offset 0. Use `run-followups.py` as the selection-policy reference; its run paths are host-specific.

Rebuild with the recorded Rust/maturin toolchains, release LTO and `extension-module,native-inchi` features. Gzipped exact dumps decompress to hashes in the source manifest. `artifact-sha256.json` indexes all other packet files. Read the decision before adopting an archived experiment.
