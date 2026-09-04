# Third-Party Notices

chematic is dual-licensed under MIT OR Apache-2.0 (see `LICENSE-MIT` /
`LICENSE-APACHE`). This file lists source-code-level derivations from
third-party projects, distinct from the compiled/linked dependency tree
covered by ordinary Cargo license metadata.

The project's copyright holder is Kentaro Tanabe (kent-tokyo). See `NOTICE` for the
project attribution notice; the third-party notices below remain separate.

## RDKit aromaticity implementation

`crates/chematic-perception/src/rdkit_parity.rs` is a source-verified,
line-level port of specific functions from RDKit's default aromaticity
model implementation, built as a diagnostic reference engine.

Derived from:

```
File:            Code/GraphMol/Aromaticity.cpp
Functions:       getAtomDonorTypeArom, countAtomElec, isAtomCandForArom,
                  applyHuckel, applyHuckelToFused
RDKit commit:    e89c9f656a694fab4105139844cba88d2e013354, an ancestor of
                 release tag Release_2026_03_4 (which resolves to
                 8afba32ec539dcb2369bc84549d802aca3f7eb39). This file is
                 byte-identical between the two commits (independently
                 diffed during Morgan M4-A0's provenance audit), so the
                 cited functions are unaffected by which commit is used.
Source URL:      https://github.com/rdkit/rdkit/blob/e89c9f656a694fab4105139844cba88d2e013354/Code/GraphMol/Aromaticity.cpp
```

```
Copyright (C) 2003-2022 Greg Landrum and other RDKit contributors
```

Licensed under the BSD 3-Clause License. Full license text below, reproduced
from RDKit's `license.txt` at the commit above.

## RDKit Morgan fingerprint hashing implementation

`crates/chematic-fp/src/rdkit_morgan_hash.rs` is a source-verified, line-level
port of RDKit's Morgan fingerprint hash-combine machinery (connectivity
invariant, bond invariant, per-round environment hash), originally built as a
diagnostic-only reference engine for Milestone M4-A0. Its core (`expand_one_pass`,
`checked_bond_invariant`) is now also reused by the production, fallible
`rdkit_morgan_ecfp4_experimental` API in `crates/chematic-fp/src/rdkit_morgan_ecfp4.rs`
(see `validation/README.md`'s Phase B section) — the diagnostic module itself
remains diagnostics-feature-gated; only the promoted internals are linked into
production.

Derived from:

```
File:            Code/GraphMol/Fingerprints/MorganGenerator.cpp
                 Code/GraphMol/Fingerprints/FingerprintUtil.cpp
                 Code/RDGeneral/hash/hash.hpp
                 Code/GraphMol/Bond.h
Functions:       MorganEnvGenerator<OutputType>::getEnvironments,
                 MorganAtomInvGenerator::getConnectivityInvariants,
                 MorganBondInvGenerator::getBondInvariants,
                 gboost::hash_combine, hash_value(std::pair), hash_range,
                 Bond::BondType enum ordinal values
RDKit release:   Release_2026_03_4
RDKit commit:    8afba32ec539dcb2369bc84549d802aca3f7eb39
Source URL:      https://github.com/rdkit/rdkit/blob/8afba32ec539dcb2369bc84549d802aca3f7eb39/Code/GraphMol/Fingerprints/MorganGenerator.cpp
```

This commit was independently resolved via the GitHub tags API
(`GET /repos/rdkit/rdkit/git/refs/tags/Release_2026_03_4`), not reused from
the aromaticity port's citation above — that citation
(`e89c9f656a694fab4105139844cba88d2e013354`) turned out to be an ancestor
commit, not the tag's actual resolution (see that section's corrected
wording above; the cited file is byte-identical between the two commits,
so no functional correction was needed there).

`crates/chematic-fp/src/morgan_environment.rs`'s own doc comment (PR #123,
predating this audit) separately cites commit
`0062b670640352ab63d6256be608615e87e1af53` for
`Code/GraphMol/Fingerprints/MorganGenerator.cpp`'s
`MorganEnvGenerator<OutputType>::getEnvironments` — independently diffed
during this audit: that commit is **not** an ancestor of the
`Release_2026_03_4` tag (diverged history, `ahead_by: 78, behind_by: 95`
from the tag resolution). The file as a whole differs by one unrelated
line (a parameter type change in a *different* function,
`updateAdditionalOutput`'s `bitId` parameter, `uint64_t`→`size_t`); the
specifically-cited `getEnvironments` function itself is byte-identical
between the two commits, so `morgan_environment.rs`'s algorithm was
unaffected — its doc comment has been updated in place with this finding
rather than left stale.

```
Copyright (C) 2003-2022 Greg Landrum and other RDKit contributors
```

Licensed under the BSD 3-Clause License (same text as above).

### Additional RDKit-derived compatibility code

The following source files also implement compatibility behavior or copied
constant/pattern data from RDKit and are covered by the same BSD-3-Clause
notice above. They must retain this notice in source and binary distributions:

```
crates/chematic-fp/src/rdkit_pattern.rs
crates/chematic-fp/src/rdkit_rdk.rs
crates/chematic-fp/src/rdkit_layered.rs
crates/chematic-fp/src/rdkit_torsion.rs
crates/chematic-fp/src/rdkit_atom_pair.rs
crates/chematic-fp/src/rdkit_isotope_delta_table.rs
crates/chematic-chem/src/gasteiger.rs
crates/chematic-perception/src/rdkit_parity.rs
```

These are compatibility implementations, not a copy or redistribution of
the RDKit build. This inventory is a copyright/license notice; it is not a
patent or freedom-to-operate opinion. No patent clearance is claimed for any
algorithm merely because its reference implementation is BSD-licensed.

## IUPAC InChI source

The optional `native-inchi` feature vendors the IUPAC InChI source under
`crates/chematic-inchi/vendor/inchi-src`. That source is distributed under
its own MIT license and retains its upstream headers. The vendored tree also
contains component notices (including the `stb_sprintf` and SHA-2 notices);
they remain applicable and are not relicensed as chematic code. This optional
feature is therefore excluded from the claim that the default build is
FFI-free.

## Research and patent boundary

Research-derived descriptors or fingerprints with unresolved FTO review are
not part of the public v1.0 surface. In particular, Spectrophores and the
rejected geometry-aware spectral fingerprint are not shipped. Existing
MMFF94, ETKDG, and other published-method implementations are independent
implementations with cited sources; this repository makes no patent-clearance
claim for them.

### BSD 3-Clause License

```
BSD 3-Clause License

Copyright (c) 2006-2015, Rational Discovery LLC, Greg Landrum, and Julie Penzotti and others
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```
