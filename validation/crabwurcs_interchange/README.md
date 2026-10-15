# crabWURCS molecular interchange check

This opt-in check covers the generic molecular boundary shared with
[crabWURCS](https://github.com/Ojas-Singh/crabWURCS): isomeric SMILES input and
MOL/SDF V2000 and V3000 output. It does not implement or validate WURCS,
glycan residue/linkage semantics, or PDB/mmCIF glycan extraction.

The Rust runner reads the fifth field of each source TSV row, parses it as
SMILES, writes both CTfile versions with a stereo-loss report, reads each block
back with diagnostics, and compares CheMatic canonical SMILES. The Python
checker uses RDKit to compare each SDF record with the corresponding source
SMILES. It classifies any difference as connectivity or stereochemistry.
Both programs exit nonzero on failure. The source tables remain in the
crabWURCS checkout and are not copied into this repository.

## Pinned source-main run, 2026-10-15

- CheMatic: `15150957a9cb4c86ffd48c0d969ff16008490602` (`origin/main` at test time).
- crabWURCS: `c83f360afc2fd38631ef3488fa7b3cce975fbba1`.
- `glycoshape_notations.tsv`: 838 rows, SHA-256 `fe9e021ac4d6788f6029d08c29170f067df5f83aba0e51da8fd63bd867327aeb`.
- `glycoshape_derived_notations.tsv`: 100 rows, SHA-256 `04710230504a03d8d47fe63facdef2fb317c648b72684cf9fccc66728bde1e34`.
- Rust: 938 source rows; zero parse, writer stereo-loss, MOL read, stereo-diagnostic, or canonical-SMILES failures for each CTfile version.
- RDKit 2025.09.3: 938/938 valid records and 938/938 matches for each CTfile version; zero connectivity or stereo differences. [Machine-readable result](rdkit-summary.json).

These are exposed source rows, not an independent holdout. Equality is limited
to the molecular graph and specified stereochemistry represented in the source
SMILES. It does not establish WURCS round trips, coverage of unseen glycans, or
behavior of a published CheMatic package.

## Reproduce

With crabWURCS checked out at the pinned commit, Rust and RDKit installed:

```sh
DATA=/path/to/crabwurcs/data
OUT=/tmp/chematic-crabwurcs-interchange
cargo run --release --manifest-path validation/crabwurcs_interchange/Cargo.toml -- "$DATA" "$OUT"
python3 validation/crabwurcs_interchange/check_rdkit.py "$DATA" "$OUT" \
  --summary "$OUT/rdkit-summary.json"
```

Check the two source hashes before comparing a rerun with the pinned result.
The SDF files are generated locally and are intentionally not committed.
