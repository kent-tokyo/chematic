# PDBQT field accuracy: bounded local comparison (2026-10-15)

Open Babel version: `Open Babel 3.2.1 -- Jul 11 2026 -- 19:21:27` (`obabel -V`).
CheMatic source: `cac480fb` plus the PDBQT parser change in this branch.
This is a two-row behavior check, not a general format-accuracy benchmark.
The exact inputs and captured Open Babel results are in the
[machine-readable record](2026-10-15-openbabel-pdbqt-fields.json).

| Input ATOM record | CheMatic `parse_pdbqt` | Open Babel `obabel -ipdbqt -ocjson` |
| --- | --- | --- |
| `ATOM      1 Cl   LIG A   1       1.000   2.000   3.000  1.00  0.00    +0.1234 Cl` | Element 17 (Cl), partial charge `0.1234` | Element 6 (C), partial charge `0.123` |
| `ATOM      1 Cl   LIG A   1       1.000   2.000   3.000  1.00  0.00           Cl` | Typed error: required partial charge is absent | Warning, then element 17 (Cl) with partial charge `0.0` |

For each row, append `\nEND\n` and pipe the text to `obabel -ipdbqt -ocjson`.
The first row deliberately uses a seven-character signed charge, wider than
the nominal six-character PDBQT charge field. Its result measures handling of
this over-width input, not strict conformance to the fixed-column specification.
The second row measures whether a missing charge is silently invented. The
CheMatic behavior is covered by `pdbqt::tests::widened_charge_preserves_precision_and_two_letter_atom_type`
and `pdbqt::tests::missing_charge_or_atom_type_is_rejected`.

Neither case establishes superiority for ordinary PDBQT records, flexible
ligand torsion trees, or other file formats. PDBQT is outside the Tier A
comparison profile in `docs/openbabel-file-io.md`.
