# Optional tetrahedral geometry check after stereo repair

`run.py` revisits the 576 exposed COD development inputs with a reported stereo
violation from `../2026-11-27-implicit-h-planarity/alternatives-after.json`.
It uses the same one deterministic ETKDG geometry, 0.8 Å new-nonbonded-pair
cutoff and at most three alternative attempts. It compares default repair
against `require_regular_geometry=True`. This is development data, not the
reserved molecular holdout. There were no processing errors, and the default
repair acceptance decisions matched all 576 previously recorded decisions.

| Default repair | With geometry check | Count |
|---|---|---:|
| Failed | Failed | 368 |
| Accepted | Rejected | 120 |
| Accepted | Accepted | 88 |
| Failed | Accepted | 0 |

Among the 88 accepted in both modes, 85 coordinate arrays were unchanged
under Python float equality. Three used alternate branches because the original
repair failed the newly required geometry check. All 88 screened outputs had
no new nonbonded pair below 0.8 Å relative to their input, and the largest
absolute bond-length change was 2.23e-15 Å. Existing close pairs may remain.

For these accepted outputs, RDKit 2026.09.2 reassigned 71 comparable declared
tetrahedral chiral tags from the same coordinates. All 71 agreed with the
screened chematic status in this fixed sample; there were no other declared
RDKit tetrahedral tags in these accepted outputs. The script verifies atom
symbols, atomic numbers, formal charges and bond-index pairs before comparing.
This checks a limited set of tetrahedral tags; it does not independently
validate E/Z stereo, unsupported stereo classes, energy, all angles or physical
conformer quality. The selected geometry screen is conservative: rejecting 120
previously accepted cases does not mean those conformers were all incorrect.

The local macOS arm64 debug wheel was built from commit `f060e04e` plus the
accompanying source changes with CPython 3.13.6 and `maturin build --locked`.
Its SHA-256 is
`42d501ac02485bc1c726185095f71a93c64747cc9cebf62ff06dd634f330ae74`. It is
not a published wheel. With that wheel and RDKit installed, reproduce using:

```sh
python validation/results/2026-11-27-stereo-repair-geometry-check/run.py > output.json
```

Results remain bounded by the verifier's stereo coverage, the selected 0.8 Å
clash cutoff, the geometry screen and the limited deterministic search. They do
not establish general RDKit parity, force-field convergence or held-out
performance.
