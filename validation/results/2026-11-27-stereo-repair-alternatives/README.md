# Bounded stereo repair alternatives on exposed development data

This experiment revisits the 576 inputs with a reported stereo violation from
`../2026-11-27-stereo-repair-clashes/results.json` using the same exposed COD
development source and one deterministic ETKDG conformer. No reserved holdout
was selected. It compares the original opt-in 0.8 Å new-pair guard against the
same guard with `max_alternative_attempts=3`. The candidate wheel is a local
macOS arm64 debug build from commit `6388b009` plus the accompanying source
changes (CPython 3.13.6, SHA-256
`f00b3f7c71033fcfa8ce9c79afcbd343e3403b2068987bf204ad6a059fd2c596`). It is
not a published package. Run with that candidate wheel and RDKit installed:

```sh
python validation/results/2026-11-27-stereo-repair-alternatives/run.py > results.json
```

There were no experiment-level exceptions. The baseline decisions reproduced
all 576 prior decisions. All 206 baseline successes remained successful with
coordinate arrays equal element by element under Python float equality. Of 370 baseline
failures, two were recovered by the bounded alternative search; the other 368
still failed. For all 208 accepted outputs, an independent Python enumeration
of nonbonded pairs found no newly introduced pair below 0.8 Å. The maximum
change in any bond length was 2.23e-15 Å. This does not mean the outputs have no
pre-existing close pairs, acceptable energies, or generally valid geometry.

RDKit 2026.09.2 reassigned tetrahedral stereochemistry from the recovered
coordinates. After checking atom symbols and bond-index pairs, both recovered
cases had one RDKit-assigned CIP center and its code agreed with the declared
input. This independent stereo check covers only those two centers.

During development, a less restrictive alternative search accepted record 719:
the internal signed-volume verifier reported success while RDKit assigned the
opposite CIP code. Its tetrahedral geometry was grossly distorted. The final
search rejects alternative geometries when a reported satisfied tetrahedral
center lies outside its four ligand positions or its ligand-direction angles
are outside 45–165°. The regression test keeps record 719 rejected. This screen
is deliberately conservative and is not an independent oracle or a guarantee
against other disagreements. The baseline accepted path is unchanged and is
not subjected to this additional alternative-only screen.

The alternatives use a shared deterministic candidate index across stereo
elements and do not enumerate all combinations. They are attempted only when
the initial repair succeeded but the close-pair guard rejected its result. The
remaining failures include both stereo-repair failures and new-pair or geometry
rejections. This small development result does not establish general conformer
accuracy, force-field convergence, or parity with RDKit.
