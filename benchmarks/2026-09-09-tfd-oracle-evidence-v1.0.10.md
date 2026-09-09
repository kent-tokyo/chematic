# RDKit TFD evidence boundary — v1.0.10

The current pinned 250-molecule comparison contains 201 finite RDKit TFD
values and 49 rigid probes with zero rotatable torsions. Those probes are
classified explicitly as `not_applicable` by the measurement harness.

Validate the checked-in evidence with:

```sh
python3 scripts/check_tfd_oracle_evidence.py
```

This closes the local TFD measurement boundary for the pinned corpus. It is
not a claim of independent force-field correctness or full conformer-quality
parity.
