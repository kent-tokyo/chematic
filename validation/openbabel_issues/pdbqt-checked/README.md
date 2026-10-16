# PDBQT fixed-width diagnostic

Original fixtures and probe authored for chematic; MIT OR Apache-2.0.
The XYZ intentionally uses coordinates too wide for PDBQT Real(8.3) fields.
No upstream attachment was copied. Open Babel #2101 reports the failure mode:
https://github.com/openbabel/openbabel/issues/2101

`openbabel_read_probe.cpp` uses the public Open Babel API with automatic bonding
disabled to report stored element/coordinate/partial-charge values independently
of chematic. It is not a docking validator. The Python runner authors accepted
and refused cases directly and compares real Python and Node/WASM output.

See `docs/pdbqt-checked-output.md` for setup, precise scope and evidence.
