# Wide Avalon fingerprints and molecular boundary contracts

This follows #813, whose processed main Codecov report is 95.58% (154,128 / 161,247 lines). The workflow command, measured files, and exclusions are unchanged.

## Avalon compatibility correction

The earlier references used 512 bits. Add 780 independent RDKit 2026.03.1 references for ten molecules at 2048, 2056, and 4096 bits, with the default flags, all flags, and each of the 24 individual flags. Include repeated nitrogen/fluorine, aromatic rings, aspirin, and linked aromatic scaffolds. All original 1,794 rows remain; the fixture now has 2,574 rows.

These exposed two implementation differences:

- Above 2048 bits, Avalon uses a seeded atom-count table. Upstream `ssmatch.c` defines `NCOUNT_SEED_HASH` as the unparenthesized `128*128`, so `seed%NCOUNT_SEED_HASH` indexes `(seed%128)*128`. The Rust port had instead used `seed%(128*128)`. Preserve the upstream slot mapping for compatibility.
- Scaffold seed products had been truncated to signed 32-bit integers before hashing. The pinned RDKit binary uses the full positive products. For benzene at 2056 requested bits (2080 internal padded bits), the scaffold-ID bit was 725 instead of 277. Compute these bounded products in 64 bits. The 2056-bit references also check padding/truncation without relying on a power-of-two modulus to hide the difference.

Sources: [RDKit Release_2026_03_1 Avalon dependency](https://github.com/rdkit/rdkit/blob/Release_2026_03_1/External/AvalonTools/CMakeLists.txt), [AvalonToolkit_2.0.5-pre.3 ssmatch.c](https://github.com/rdkit/ava-formake/blob/AvalonToolkit_2.0.5-pre.3/src/main/C/common/ssmatch.c). Regenerate with `scripts/generate_avalon_boundary_fixtures.py` and exactly RDKit 2026.03.1. The reference comparison establishes parity with that pinned binary, not with arbitrary compiler behavior for upstream signed-overflow expressions.

## Additional contracts

- InChI reading: extract all 37 unique InChI literals from RDKit's fixed-version official `External/INCHI-API/test.cpp`, including adjacent C++ string fragments. Compare canonical output or rejection with `Chem.MolFromInchi`. Retain unusual valence, isotopes, radicals, stereochemistry, and the three rejected inputs. `scripts/generate_inchi_reader_boundary_fixtures.py` verifies the upstream source SHA-256 (`d3faba4a81f82c03647c275f0c0978cd7887e883b0fd9a3ca3e473252d1a9b01`) before generating the TSV; pass that source file as its argument.
- Semantic expansion: exact atom budgets, source-to-expanded atom identities, one through three repeat units, and both legacy and typed end caps. Compare expanded graphs with explicit cyclic amine/alcohol structures. Reject malformed endpoints, wildcard topology, end caps, repeat metadata, and R-group alternatives without editing the source model or molecule.
- Stereo annotations: contradictory wedge directions identify the original centre; redundant CF4 and undercoordinated tagged carbon produce distinct diagnostics, retaining their original tags and bonds.
- Parent selection: zero-timeout results retain the input and an empty transform audit; a zero transform budget reports unfinished nitroso-to-oxime conversion. Cytosine and guanine keto/enol forms converge to idempotent parents with bounded audit atom/bond identities. These cases do not claim full tautomer enumeration.
- Molecular alignment: permuted equivalent ammine ligands around Cu map without reversing dative bonds. Pinned RDKit gives CalcRMS/GetBestRMS zero and the default first-map fit 0.9307312929830155. Reject invalid atom maps and weight counts. Check mapped residuals and Cu identity.
- Sphere overlap: independent equal-unit-sphere lens volume 5π/12 at separation 1, full containment 4π/3, tangency/disjoint zero, exchange symmetry, and cubic scaling. Inputs have positive centre separation; this does not establish the coincident-centre helper behavior.
- Legacy torsion API: enumerate actual four-atom paths on 23 molecules, with implicit/explicit hydrogens. Preferences produce finite nonnegative, 360-degree-periodic penalties and zero at the stated angle; both preference/None paths occur. This checks the retained numerical API, not the physical validity of its heuristic angles; see `docs/rfcs/3d_torsion_knowledge_audit.md`.
- Reaction export: a valid rich document reports combined step/provenance/condition/coefficient loss categories without changing the document, while plain reaction SMILES round-trips.
- SVG grid: highlights stay in their cell, invalid atom indices are ignored, per-atom/default colors are escaped and sorted deterministically, highlighted bonds stay in the first cell, and coordinates remain finite.

Coverage is execution evidence for these specific contracts, not a guarantee of general chemical correctness. Validation results and the processed Linux service report are recorded on the accompanying PR.
