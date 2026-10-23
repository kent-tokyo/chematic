# Contributor corpus policy

Every checked-in molecule, reaction, expected-value fixture, benchmark input,
or compatibility corpus needs enough provenance to be redistributed and
regenerated safely.

## Required metadata

Record a stable local ID, original source and retrieval date, license or
`NO-REDISTRIBUTION` decision, redistribution status, minimization rationale,
pinned oracle and options, expected result with units/tolerance, and SHA-256 of
the normalized payload.

For `derived-only` data, commit the recipe and hash rather than the source
payload. For excluded material, use a synthetic or minimized replacement.
Personal, patient, proprietary, or credential-bearing data must not enter the
repository.

## Regeneration and review

- Pin oracle versions and distinguish defaults from explicit options.
- Never overwrite expected values silently after an oracle upgrade; record the
  new oracle and explain the difference.
- Normalize encoding, line endings, and numeric formatting before hashing.
- A clean checkout must be able to regenerate the record from documented,
  non-secret inputs.
- Review provenance, licensing, minimization, and unnecessary third-party data
  separately from test correctness.

A passing test proves only the checked-in assertion. It does not establish
license clearance or scientific validity beyond the recorded oracle contract.
