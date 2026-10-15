# Default TPSA zero-contribution guard evidence

Read decision.json, measurement-plan.json, source-manifest.json, both frontier matrices and adaptive-selection.json. Any focused matrices are separate follow-ups. All blocks, modes, gates and hashes are retained.

Extract candidate-build-source.tar.gz to obtain source/. Replace paths with baseline-overrides/ and remove the added test listed in the manifest to reconstruct the immediate baseline. Source-reconstruction.json checks every hash. Read the adoption decision before using an archived experiment.

The wheel directories contain exact measured binaries. Install separately with Python 3.13 and pip install --no-deps --target <directory> <wheel>. Comparator dependencies are pinned in comparator-environment.txt. Build with recorded release LTO, extension-module/native-inchi features, Rust and maturin versions.

Run the extracted scripts/bench_python_speed_frontier.py using --corpus source/validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi --offset 5000 --limit 1000 --blocks 20, and supply baseline/candidate wheel and installation paths, source manifest, base revision and output path. Repeat at offset 0. run-frontiers.py and run-followups.py record exact host commands and selection policy; paths need adaptation on another host. Gzipped exact dumps match manifest hashes. artifact-sha256.json indexes all other packet files.
