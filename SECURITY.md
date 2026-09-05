# Security Policy

## Supported versions

| Version | Supported | Status |
|---|---|---|
| v1.0.7 | Yes | Current release |
| v1.0.6 | Security fixes only | Upgrade to v1.0.7 is recommended |
| v1.0.0-v1.0.5 | Security fixes only | Upgrade to v1.0.7 is recommended |
| v0.89.0 | Security fixes only | Previous published release |
| Earlier versions | No | End of life |

The current release is published on GitHub Releases, crates.io, PyPI, and npm.
Release history belongs in [`CHANGELOG.md`](CHANGELOG.md); detailed historical
security notes are archived in
[`docs/archive/security-policy-through-v1.0.4.md`](docs/archive/security-policy-through-v1.0.4.md).

## Reporting a vulnerability

Use [GitHub private vulnerability reporting](https://github.com/kent-tokyo/chematic/security/advisories/new)
for non-public reports. For sensitive coordination, contact
`36805997+kent-tokyo@users.noreply.github.com`.

Please include the affected package/version, impact, minimal reproducer, and
whether untrusted input or an optional feature is required. Do not open a
public issue before coordinated disclosure when exploitation details are
sensitive.

Expected handling:

- initial response within seven days;
- target fix within 30 days, depending on severity and reproducibility;
- coordinated disclosure after a patched artifact is available;
- reporter credit unless anonymity is requested.

## Security boundary

In scope:

- memory unsafety or undefined behavior;
- unbounded CPU, memory, recursion, or output growth on public untrusted-input
  paths;
- parser or serialization behavior that crosses a trust boundary;
- command, path, network, credential, or supply-chain injection;
- bypass of documented WASM, MCP, CLI, or binding resource limits.

Usually handled as correctness issues rather than vulnerabilities:

- chemistry accuracy or compatibility differences without a security impact;
- unsupported RDKit/Open Babel behavior;
- Experimental 3D/MMFF94 quality or convergence limitations;
- model, method, or dataset quality claims.

If a correctness defect can cause unsafe downstream action, denial of service,
or a trust-boundary violation, report it privately so impact can be assessed.

## Runtime and dependency model

- The common Rust, Python, and WASM chemistry paths are pure Rust and perform
  no implicit network access.
- `chematic-mcp` is mostly local, but its explicit `name_to_smiles` tool calls
  the PubChem REST API. Deployments must apply their own egress, privacy,
  timeout, and availability policy to that tool.
- The optional `native-inchi` feature uses a reviewed vendored C FFI boundary
  and is unavailable to WASM. The default pure-Rust InChI path is approximate.
- Filesystem access occurs only through APIs or commands that the caller
  explicitly invokes. Applications remain responsible for path permissions,
  quotas, and sandboxing.
- Browser builds run inside the browser sandbox. Three-browser smoke and
  adversarial tests cover representative malformed and oversized inputs, but
  do not constitute a complete browser security audit.

## Maintainer controls

Repository-local controls include:

- finite parser, search, binding, MCP, CLI, and serialization limits with
  typed failures;
- four fuzz targets and minimized regression corpora;
- focused Miri and Linux ASan/LSan/TSan lanes;
- Rust/Python/Node/WASM contract tests and Chromium/Firefox/WebKit smoke tests;
- dependency and license checks, immutable GitHub Action pins, SBOM,
  provenance, checksums, and release-key verification;
- an unsafe-surface allowlist isolating the optional native-InChI boundary.

Reproduction commands and evidence boundaries are documented in
[`docs/v1.0-local-release-gate.md`](docs/v1.0-local-release-gate.md),
[`docs/security-surface.md`](docs/security-surface.md), and
[`docs/security-review/`](docs/security-review/). Independent non-maintainer
review remains an external follow-up and is not claimed as completed.

## User guidance

- Use the latest patch release and pin dependencies and artifacts where
  reproducibility matters.
- Treat all molecular and document input as untrusted; keep default limits or
  choose explicit lower limits for exposed services.
- Validate chemistry output for the intended scientific or regulatory use.
- Keep Experimental 3D/MMFF94 output behind the documented sanity and failure
  checks.
- Restrict or disable the PubChem-backed MCP tool when network access or data
  disclosure is not acceptable.

Copyright and third-party notices are in [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md).
Implementation and patent/FTO boundaries are in
[`docs/implementation-provenance.md`](docs/implementation-provenance.md).

Last updated: 2026-09-05.
