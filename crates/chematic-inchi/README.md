# chematic-inchi

InChI and InChIKey generation for Rust. Two modes:

| Mode | Function | Compliance | WASM |
|------|----------|-----------|------|
| **Standard** (`native-inchi` feature) | `standard_inchi` / `standard_inchi_key` | Bit-exact with IUPAC InChI 1.07.5 | No (requires C compiler) |
| **Approximate** (default) | `inchi` / `inchi_key` | Non-standard — different canonical order, no tautomer normalization | Yes (zero C deps) |

## Standard IUPAC InChI (`native-inchi` feature)

Enable the `native-inchi` feature to link the vendored IUPAC InChI C library (v1.07.5, MIT).
Output is **bit-exact** with the reference implementation.

```toml
[dependencies]
chematic-inchi = { version = "1.0.0", features = ["native-inchi"] }
```

```rust
use chematic_inchi::{standard_inchi, standard_inchi_key};
use chematic_smiles::parse;

let mol = parse("c1ccccc1").expect("benzene");
let inchi_str = standard_inchi(&mol).unwrap();
// Produces the IUPAC-standard InChI (canonical order matches reference):
// InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H
let key = standard_inchi_key(&inchi_str).unwrap();
assert_eq!(key, "UHOVQNZJYSORNB-UHFFFAOYSA-N");
```

All layers are handled by the IUPAC reference library: formula, connectivity, hydrogen,
tautomer normalization, mobile-H, stereo (/b E/Z double-bond, /t tetrahedral, /m, /s),
charge (/q, /p), and isotope (/i).

> **Known limitation:** 2 of 5,000 drug-like molecules (0.04%) still fail with
> `InchiError::KekulizationFailed`: one boron-aromatic ring (`b1ccccn1`) and pure H₂
> (`[H][H]`, which has no heavy atoms).  All other previously failing molecules —
> including indolizine-type bridgehead-N systems (~109 cases) and complex polycyclics
> (~17 cases) — are now handled by the 4-pass kekulization algorithm in `chematic-core`.
> No silent wrong output is produced; the error is always explicit.

## Approximate InChI (WASM / zero-C-dependency fallback)

The default `inchi()` / `inchi_key()` functions are **not standard IUPAC InChI**. They use
a different canonical ordering and do not implement tautomer or mobile-H normalization.
Use them only when WASM compatibility or a zero-C-dependency build is required.

```rust
use chematic_inchi::inchi;
use chematic_smiles::parse;

let mol = parse("c1ccccc1").expect("benzene");
// Non-standard approximation — does NOT match the IUPAC reference:
let approx = inchi(&mol);
```

## Dependencies

- [`chematic-core`](../chematic-core/README.md) — molecular graph types
- [`chematic-smiles`](../chematic-smiles/README.md) — SMILES parser
- [`chematic-chem`](../chematic-chem/README.md) — CIP rules for stereo assignment
- [`sha2`](https://docs.rs/sha2/) — SHA-256 for approximate InChIKey hash
- `cc` *(optional, `native-inchi` only)* — C compiler for vendored IUPAC library

## References

- IUPAC InChI standard: https://iupac.org/what-we-do/inchi/
- InChI Trust: https://www.inchi-trust.org/

## See Also

- [`chematic`](../chematic) — main umbrella crate
- [`chematic-smiles`](../chematic-smiles) — SMILES I/O
- [`chematic-depict`](../chematic-depict) — SVG 2D structure drawing
