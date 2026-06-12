# Security Policy

## Supported Versions

| Version | Supported | Status |
|---------|-----------|--------|
| v0.1.89 | Yes | Current release — active support |
| v0.1.88 | Limited | Security fixes only (limited) |
| < v0.1.87 | No | Unsupported |

**Active Support**: Latest release (v0.1.89) receives all security updates.  
**Limited Support**: Previous release receives critical security fixes only.  
**End of Life**: Older versions receive no support.

---

## Reporting a Vulnerability

### For Non-Public Reports (Recommended)

Please use **GitHub Security Advisories** to report vulnerabilities privately:

1. Go to: https://github.com/kent-tokyo/chematic/security/advisories
2. Click **Report a vulnerability** (or **New draft security advisory**)
3. Fill in the vulnerability details:
   - **Vulnerability type**: Choose from the dropdown
   - **Package**: chematic (or specific sub-crate)
   - **Severity**: High / Medium / Low
   - **Description**: Clear reproduction steps and impact
4. Submit the draft

**GitHub will**:
- Notify maintainers immediately
- Allow collaborative discussion in a private space
- Coordinate embargo if needed before public disclosure

### For Sensitive Matters

Email: **36805997+kent-tokyo@users.noreply.github.com**

---

## Response Timeline

- **Initial Response**: Within 7 days of report submission
- **Security Fix**: Aim for 30 days or less
- **Disclosure**: Coordinated disclosure after patch release (usually immediate)
- **Credit**: Vulnerability reporter will be credited unless they request anonymity

---

## Scope

### In Scope

Vulnerabilities affecting the security of chematic itself:

- **Supply chain attacks**: Malicious dependencies (Dependabot monitors these)
- **Memory safety**: Unsafe code blocks that could lead to UB
- **Cryptographic flaws**: In hashing or data serialization
- **Input validation**: Buffer overflows, panic on invalid SMILES/InChI
- **Privilege escalation**: Not applicable to a library, but relevant for WASM sandbox

### Out of Scope

The following are **NOT** considered security vulnerabilities in chematic:

- **Incorrect chemistry results**: Accuracy bugs (report as GitHub Issues instead)
- **Missing RDKit features**: Partial implementation vs. RDKit (by design)
- **Dependency vulnerabilities**: If a transitive dependency has a CVE, report to that project; we will update via Dependabot
- **Third-party FFI exploits**: chematic has zero C/C++ FFI by design
- **Transition metal chemistry**: Out of scope (atom valence model limitation, not a vulnerability)
- **ML model attacks**: Cheminformatics is non-ML in chematic

---

## Security Best Practices for Users

### Dependency Updates

This project uses **Dependabot** to automatically track security updates:

- Cargo dependencies: Updated weekly
- GitHub Actions: Updated weekly
- All updates create pull requests for review before merge

Subscribe to GitHub notifications for this repository to receive alerts about security updates.

### Using chematic Safely

1. **Keep chematic updated**: Use latest published version from npm/crates.io
2. **Validate chemistry output**: This is a chemistry library — always validate results are chemically sensible
3. **No network access**: chematic has zero network calls; it's safe for sandboxed/offline environments
4. **No file I/O side effects**: WASM and Rust versions are side-effect-free (except explicit file operations you request)

---

## Security Considerations by Use Case

### Browser (WASM)

- Safe: Runs in browser sandbox
- Safe: No network calls
- Note: SMILES/InChI parsing is complex; malformed input won't exploit chematic but may consume CPU

### Node.js / Electron

- Safe: Same as browser, plus npm registry security
- Dependabot monitors npm dependencies
- Note: Electron-specific sandbox rules apply (not chematic's responsibility)

### Rust / Server

- Safe: No FFI, no network calls, no file I/O unless you explicitly request it
- Dependabot monitors cargo.io registry
- Note: If you use unsafe code interop with chematic, validate chemical outputs before use

---

## Acknowledgments

Security researchers who responsibly disclose vulnerabilities help keep chematic safe for everyone. We appreciate their efforts and will credit them publicly (unless they prefer anonymity).

---

## GitHub Security Settings

This repository has the following GitHub security features enabled:

- Secret scanning: Alerts if credentials are accidentally pushed
- Dependabot alerts: Weekly dependency updates
- Dependabot security updates: Automatic PR creation for vulnerabilities
- Code scanning: Cargo audit in CI (via GitHub Actions)
- Private vulnerability reporting: This security policy file enables it

---

**Last Updated**: 2026-06-12  
**Security Contact**: 36805997+kent-tokyo@users.noreply.github.com
