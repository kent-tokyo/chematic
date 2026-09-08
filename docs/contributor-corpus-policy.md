# Contributor corpus policy

This policy applies to every checked-in molecule, reaction, expected-value
fixture, benchmark input, and compatibility corpus. The goal is a small,
reproducible corpus that can be redistributed and regenerated without silently
turning a third-party dataset into project data.

## Required record

Each new corpus record or manifest entry must identify:

| Field | Requirement |
|---|---|
| `id` | Stable local identifier; do not encode a temporary filename. |
| `provenance` | Original source, URL or citation, retrieval date, and transformation steps. |
| `license` | SPDX identifier or an explicit `NO-REDISTRIBUTION` decision. |
| `redistribution` | `allowed`, `derived-only`, or `excluded`, with a short rationale. |
| `minimization` | Why the smallest structure/input sufficient for the assertion was selected. |
| `oracle` | Engine/tool name, exact version, and invocation or settings. |
| `expected` | Expected value/status plus units and tolerance where applicable. |
| `sha256` | Hash of the normalized fixture payload, not a mutable source archive. |

For `derived-only` data, check in the recipe and hash, not the source payload.
For `NO-REDISTRIBUTION` or `excluded` material, keep only a synthetic or
minimal replacement that tests the same behavior. Personal, patient,
proprietary, or credential-bearing data must never enter the repository.

## Oracle and regeneration rules

- Pin oracle versions in the manifest and in the benchmark/validation note.
- Record whether the oracle was run with defaults or explicit settings.
- Do not silently refresh expected values after an oracle upgrade. Add a new
  oracle version, explain the delta, and retain the old result when it remains
  useful for compatibility.
- Normalize line endings, record encoding, and numeric formatting before hashing.
- A reviewer must be able to regenerate a record from a clean checkout using
  documented, non-secret inputs. Network-only regeneration is not a release
  gate unless the source and a cached hash are both documented.

## Review gate

The contributor and reviewer must confirm provenance, licensing, minimization,
oracle pinning, and that the fixture contains no unnecessary identifying or
third-party bulk data. CI may reject a record missing any required field. A
passing local test proves only the checked-in assertion; it does not establish
license clearance or scientific validity beyond the pinned oracle contract.
