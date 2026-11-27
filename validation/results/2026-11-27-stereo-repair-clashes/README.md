# Stereo repair: pairwise clash rejection on exposed development data

`run.py` reads the 1,345-row exposed COD development input from
`../2026-11-27-selenium-crystal-stereo/source-development-input.json` (SHA-256
`e3c2cf4aa5c7c5740fd3d23c670e26535515b97f16a02ee2bf387661513726a1`). This is
development data, not a reserved holdout. It skips strings without `@`, `/` or
`\\`; marker presence does not imply the verifier can evaluate all specified
stereo. For the remainder it generates one deterministic ETKDG coordinate set,
selects cases with a reported violation, and compares plain repair with the
optional `clash_cutoff=0.8` repair. It independently enumerates pair distances
in Python using the same molecule's bond table and compares the sets of new
pairs and the guard decision. No experiment-level exceptions occurred.

The local macOS arm64 debug wheel was built from commit `afe8194e` plus the
accompanying source changes, using CPython 3.13.6 and `maturin build --locked`.
Its SHA-256 is
`8631ae32b8a210a3bcb9e1bd97facfa7bcb2b0a1f0173856993e981ece406fce`. It is not a
published wheel. Reproduce with the candidate wheel installed:

```sh
python validation/results/2026-11-27-stereo-repair-clashes/run.py > results.json
```

Of the 1,345 rows, 539 had no selected marker and 230 had no reported stereo
violation. Among the remaining 576, ordinary stereo repair failed for 256;
206 succeeded and passed the optional guard, while 114 succeeded without the
guard but were rejected because they introduced one or more new close pairs.
The independent pair enumeration and guard agreed for all 320 successful
ordinary repairs. In 14 rejected cases the total number of close pairs stayed
the same or decreased, illustrating why counts alone are insufficient.

The screen uses an absolute 0.8 Å cutoff, includes 1–3 pairs and separate
fragments, and only rejects pairs that newly cross below the cutoff. It does not
reject a pre-existing close pair merely because it gets closer. The 206 accepted
geometries may already contain clashes; this evaluation does not establish
physical quality, force-field convergence, general accuracy, or RDKit parity.
