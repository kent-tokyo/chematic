# Neutral selenium stereo: fixed-coordinate development comparison

This probe evaluates a limited class of neutral three-coordinate selenium:
no implicit hydrogen, two ordinary single bonds, and exactly one double bond
to N, O, S or Se. The new verifier uses the existing lone-pair direction
convention for this class. Other valence, charge and partner cases are not
claimed as supported.

Six authored development SMILES patterns use both `@` and `@@`, seeds 7 and
42, and multiple atom orders. Six additional authored control patterns use
both tags and seeds 13 and 99. All are small synthetic cases chosen during
this work; neither group is an independent holdout or a representative
population. RDKit 2026.09.2 generated coordinates with explicit hydrogens.
For each generated geometry and its reflected image, the probe checks that
the indexed atom/bond graph matches, then independently calls RDKit's 3D
stereo reassignment and compares its resulting selenium CIP label with the
input declaration. It also calls chematic's boolean stereo gate on the exact
same coordinates.

| Group | Coordinates, including mirrors | RDKit agrees with declaration | Changed chematic gate agrees | Before gate returned true |
| --- | ---: | ---: | ---: | ---: |
| Development | 48 | 24 | 48/48 | 0/48 |
| Authored controls | 48 | 24 | 48/48 | 0/48 |

The old boolean gate returned false for every case; it does not distinguish
an unevaluable center from a violated one. The changed gate returned true for
all 48 original geometries and false for all 48 mirrors, matching RDKit's 3D
reassignment in this sample. The coordinates matched between the before and
after probe runs. The Rust regression fixture separately asserts satisfied,
violated and unevaluable statuses, including degenerate and unsupported
coordination; its 24 forward rows are the exposed development examples.

The previously exposed square-planar selenium case also changed from
`other_stereo` unevaluable to a returned checked embedding at seed 42 and
planarity threshold 0.3. RDKit's 3D reassignment reported the same `S` CIP
label as the input for the other selenium center; its reflected coordinates
were independently assigned `R`. This one case does not establish the
correctness of every square-planar assignment or any bond, energy, physical
stability, or crystal-structure claim.

RDKit emitted UFF typing warnings for some selenium examples; the scripted
run suppresses those warnings. Coordinates are used here only as parity test
inputs, not as validated physical geometries. The controls were authored
alongside the change and do not prove general accuracy.

The host used macOS 27.0.1 arm64, Python 3.13.6, RDKit 2026.09.2 and
`RAYON_NUM_THREADS=1`. Both CPython wheels were nonpublished source builds
with release LTO disabled and 16 codegen units. Before wheel SHA-256:
`65740f3f6d3d596c26001977d55f0078f98135f88c70b70979702ae75cd0ac1d`;
after wheel SHA-256:
`87abca4c0637601200ec01fbbc62e7846a137227d0f096eee348b2832f12ea3c`.
The after extension SHA-256 is
`63a9669b030998b6f58ab930426648423979d905e0364ac489c17e8897789259`.
The [comparison script](../scripts/validation/neutral_selenium_stereo_development.py)
SHA-256 is `86e1149e39f074bb42ad8c298920ecb0281522ae72d21db888cdd3d2c5c9a83b`.
The [checked square-planar probe](../scripts/validation/neutral_selenium_checked_square_planar.py)
SHA-256 is `a4eff5be5a5fdc42badbba87dcc96f46b484d5103482c4850bf4803056e90439`.
The [fixed square-planar oracle generator](../scripts/validation/neutral_selenium_square_planar_oracle.py)
SHA-256 is `4f6c7dfb47df92b5043bee7472ba1b1aba57f884f2f6b75ba6d600b01d2b8016`.
Raw inputs, coordinates, CIP labels and gate results are retained in
`validation/results/2026-11-27-neutral-selenium-stereo/`.
