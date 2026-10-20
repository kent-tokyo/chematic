# MMFF dielectric update contract

`cases.json` is original fixture data: RDKit 2026.09.1, explicit-H SMILES,
`EmbedMolecule(randomSeed=42)` for organic molecules, manually placed ions,
constant dielectric model 1, 100 Å nonbonded cutoff and interfragment terms
enabled. It contains 56 energy/gradient references for seven cases, two MMFF
variants and epsilon 1, 4, 78 and 0.5. No Open Babel attachment is copied.

The original C++ probe reproduces Open Babel #2811 using Na⁺/Cl⁻ at 4 Å.
For a Homebrew installation:

```sh
c++ -std=c++17 openbabel_probe.cpp -I/opt/homebrew/include/openbabel3 -L/opt/homebrew/lib -lopenbabel -o /tmp/openbabel-dielectric-probe
```

From the repository root, with the newly built chematic Python module and
RDKit 2026.09.1 on `PYTHONPATH`:

```sh
python3 scripts/validation/validate_mmff_dielectric.py --openbabel-probe /tmp/openbabel-dielectric-probe --output /tmp/mmff-dielectric-evidence
```

The raw packet separates the stale-cache Open Babel diagnostic from the
native-versus-RDKit energy/gradient/minimization comparisons.
