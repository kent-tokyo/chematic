# Third-Party Notices

chematic is dual-licensed under MIT OR Apache-2.0 (see `LICENSE-MIT` /
`LICENSE-APACHE`). This file lists source-code-level derivations from
third-party projects, distinct from the compiled/linked dependency tree
covered by ordinary Cargo license metadata.

## RDKit aromaticity implementation

`crates/chematic-perception/src/rdkit_parity.rs` is a source-verified,
line-level port of specific functions from RDKit's default aromaticity
model implementation, built as a diagnostic reference engine (see
`docs/aromaticity_a1_rfc.md`).

Derived from:

```
File:            Code/GraphMol/Aromaticity.cpp
Functions:       getAtomDonorTypeArom, countAtomElec, isAtomCandForArom,
                  applyHuckel, applyHuckelToFused
RDKit release:   Release_2026_03_4
RDKit commit:    e89c9f656a694fab4105139844cba88d2e013354
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
invariant, bond invariant, per-round environment hash), built as a
diagnostic-only reference engine for Milestone M4-A0 (see
`docs/` and the PR description for the design writeup). Not wired into any
production API.

Derived from:

```
File:            Code/GraphMol/Fingerprints/MorganGenerator.cpp
                 Code/GraphMol/Fingerprints/FingerprintUtil.cpp
                 Code/RDGeneral/hash/hash.hpp
Functions:       MorganEnvGenerator<OutputType>::getEnvironments,
                 MorganAtomInvGenerator::getConnectivityInvariants,
                 MorganBondInvGenerator::getBondInvariants,
                 gboost::hash_combine, hash_value(std::pair), hash_range
RDKit release:   Release_2026_03_4
RDKit commit:    8afba32ec539dcb2369bc84549d802aca3f7eb39
Source URL:      https://github.com/rdkit/rdkit/blob/8afba32ec539dcb2369bc84549d802aca3f7eb39/Code/GraphMol/Fingerprints/MorganGenerator.cpp
```

This commit was independently resolved via the GitHub tags API
(`GET /repos/rdkit/rdkit/git/refs/tags/Release_2026_03_4`), not reused from
the aromaticity port's citation above — that citation
(`e89c9f656a694fab4105139844cba88d2e013354`) turned out to be an ancestor
commit, not the tag's actual resolution; not corrected here, flagged for a
future cleanup.

```
Copyright (C) 2003-2022 Greg Landrum and other RDKit contributors
```

Licensed under the BSD 3-Clause License (same text as above).

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
