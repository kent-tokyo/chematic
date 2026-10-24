# CML lexical meaning packet

Original synthetic MIT fixtures authored on 2026-10-24: three controls and
thirteen reproductions of graph corruption caused by XML layout or opaque
content. `cases.json` is the authored graph/H/2D-coordinate oracle.

`baseline.json` and `candidate.json` retain every route's output or refusal,
runtime identities and observations. Their graph/H comparison is separate from
source-byte retention and native coordinate tests. Baseline source commit:
`342d1d5f9ad09aa6ade4e3079a902cba856209e2`; affected CML sources match main base
`fcbf6d0a`. Candidate source hashes identify the uncommitted source snapshot.

See [scope, counts and reproduction](../../docs/reference/cml-lexical-meaning.md).
This development packet reaches Open Babel 3.2.1's result; it does not establish
general CML, XML schema, stereochemistry or speed superiority.
