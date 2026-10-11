# Indigo workflow comparison

This lane compares the bounded CheMatic reaction, standardization, and R-group
surfaces with Indigo's public Python API. It is a development gate, not a
toolkit-wide superiority claim.

The runner uses:

- the checked-in 83-row stratified reaction corpus, with RDKit 2026.03.6 as
  the product-set reference;
- the checked-in Phase 1 standardization fixtures and their declared expected
  largest fragments;
- a checked-in labelled terminal-wildcard R-group matrix with explicit
  substituent payloads.

Timing is reported only for operations where both engines satisfy the declared
accuracy observer. Twenty-one alternating-order blocks are used by default.

```bash
.venv/bin/python validation/indigo_workflows/benchmark.py \
  --output benchmarks/indigo-workflows-source.json
```

The required comparator is `epam.indigo==1.46.0`. The runner records exact
engine versions, source commit, corpus hashes, row accounting, errors, and raw
timing blocks.
