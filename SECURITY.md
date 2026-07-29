# Security Policy

## Supported Versions

| Version | Supported | Status |
|---------|-----------|--------|
| v0.8.0 | Yes | Current release — active support |
| v0.4.22 | Limited | Security fixes only (limited) |
| < v0.4.22 | No | Unsupported |

**Active Support**: Latest release (v0.8.0) receives all security updates.  
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

## Security Fix History

### v0.4.14 (2026-06-21)

Seven security-relevant fixes shipped in this release:

1. **KET parser: R-element silently mapped to carbon** — `starts_with('R')` was matching real elements (Ru, Rh, Re, Rn) and silently treating them as carbon. Fixed to only match actual R-group labels (R, R#, R1–R9), ensuring legitimate elements are parsed correctly.

2. **KET parser: no atom/bond count limit (Memory-DoS)** — A crafted large KET file could allocate unbounded memory. Fixed with hard guards: `MAX_ATOMS = 10,000` and `MAX_BONDS = 20,000`; inputs exceeding these limits are rejected with an error.

3. **KET parser: isotope u64→u16 silent truncation** — Isotope values larger than 65,535 were silently truncated to incorrect values. Fixed with an explicit bounds check that returns a parse error for out-of-range isotopes.

4. **RInChI: malformed separator for empty products** — When a reaction had no products, the `!d-` separator token appeared as a second reactant, producing an incorrect InChI string. Fixed by only emitting the product block when it is non-empty.

5. **SMIRKS transform: stereo_groups dropped** — `clear_orphaned_stereo_bonds` discarded ABS/OR/AND stereo groups on the transformed molecule. Fixed by restoring stereo groups via `set_stereo_groups()` after the orphan-bond sweep.

6. **WASM docs: wrong GETAWAY vector size** — Public API documentation claimed 9 values for the GETAWAY descriptor (actual: 19) and 19 for the combined vector (actual: 29). Corrected docstrings prevent callers from allocating undersized buffers.

7. **Hosoya index: exponential recursion (CPU-DoS)** — No guard on molecule size allowed the public API to hang on large inputs due to exponential recursive matching. Fixed with an `n > 40` sentinel that returns early, bounding computation.

### Prior Releases

- **Prng fixed seed (Weyl counter)** — PRNG seeding used a fixed Weyl counter, reducing entropy. Fixed to use a properly varied seed source.
- **MCP DoS** — `find_mcs` lacked a timeout and atom-count limit, and `name_to_smiles` did not URL-encode user input. Both fixed with appropriate guards and encoding.
- **CIF parser safety** — Division-by-zero and missing unit-cell parameter checks in the CIF parser could panic on malformed crystallographic files. Fixed with explicit validation.
- **SMARTS atom map fixes** — Atom-map handling in SMARTS patterns produced incorrect match results in edge cases. Fixed to correctly propagate atom maps through the VF2 isomorphism engine.

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

**Last Updated**: 2026-06-27  
**Security Contact**: 36805997+kent-tokyo@users.noreply.github.com
