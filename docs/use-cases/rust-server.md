# Use case: Chemistry REST API in Rust (Axum)

## Problem

You're building a backend service in Rust and need chemistry capabilities — descriptor calculation, similarity search, substructure filtering. Calling out to a Python subprocess or embedding RDKit adds a C++ build dependency and a foreign runtime to an otherwise single-binary deployment.

## Solution

chematic is a pure-Rust library. `cargo add chematic` is all you need — no Python runtime, no shared C++ libraries, no Emscripten. The same binary that serves HTTP also does chemistry.

## Output / What you get

```bash
$ curl -s -X POST http://localhost:3000/descriptors \
    -H "Content-Type: application/json" \
    -d '{"smiles": "CC(=O)Oc1ccccc1C(=O)O"}'

{"mw":180.16,"logp":1.31,"tpsa":63.6,"hbd":1,"hba":4,"lipinski_passes":true}
```

## Setup

```toml
# Cargo.toml
[dependencies]
chematic = { version = "1.0.4", features = ["chem", "fp", "smarts"] }
axum      = "0.7"
serde     = { version = "1", features = ["derive"] }
tokio     = { version = "1", features = ["full"] }
```

## Descriptor endpoint

```rust
use axum::{extract::Json, routing::post, Router};
use chematic_smiles::parse_smiles;
use chematic_chem::{molecular_weight, logp_and_mr, tpsa, hbd_count, ring_bundle};
use serde::{Deserialize, Serialize};

#[derive(Deserialize)]
struct SmilesInput { smiles: String }

#[derive(Serialize)]
struct DescriptorOutput {
    mw: f64, logp: f64, tpsa: f64, hbd: usize, hba: usize,
    lipinski_passes: bool,
}

async fn descriptors(Json(body): Json<SmilesInput>) -> Json<DescriptorOutput> {
    let mol = parse_smiles(&body.smiles).expect("invalid SMILES");
    let rb = ring_bundle(&mol);
    let mw = molecular_weight(&mol);
    let (logp, _) = logp_and_mr(&mol);
    let tpsa = tpsa(&mol);
    let hbd = hbd_count(&mol);
    Json(DescriptorOutput {
        mw, logp, tpsa, hbd, hba: rb.hba_count,
        lipinski_passes: mw <= 500.0 && hbd <= 5 && rb.hba_count <= 10 && logp <= 5.0,
    })
}

#[tokio::main]
async fn main() {
    let app = Router::new().route("/descriptors", post(descriptors));
    axum::serve(
        tokio::net::TcpListener::bind("0.0.0.0:3000").await.unwrap(),
        app,
    ).await.unwrap();
}
```

```bash
curl -X POST http://localhost:3000/descriptors \
  -H "Content-Type: application/json" \
  -d '{"smiles": "CC(=O)Oc1ccccc1C(=O)O"}'

# {"mw":180.16,"logp":1.31,"tpsa":63.6,"hbd":1,"hba":4,"lipinski_passes":true}
```

## Similarity search endpoint

```rust
use chematic_fp::ecfp4;
use chematic_smiles::parse_smiles;

async fn similarity(Json(body): Json<SimilarityInput>) -> Json<Vec<Hit>> {
    let query_mol = parse_smiles(&body.query).unwrap();
    let query_fp  = ecfp4(&query_mol);

    let hits: Vec<Hit> = body.library
        .iter()
        .enumerate()
        .filter_map(|(i, smi)| {
            let mol = parse_smiles(smi)?;
            let fp  = ecfp4(&mol);
            let sim = chematic_fp::tanimoto_bitvec(&query_fp, &fp);
            (sim >= body.threshold).then(|| Hit { index: i, smiles: smi.clone(), score: sim })
        })
        .collect();

    Json(hits)
}
```

## SDF batch processing

```rust
use chematic_mol::SdfReader;
use std::io::BufReader;

fn process_sdf(path: &str) -> Vec<DescriptorOutput> {
    let f = std::fs::File::open(path).unwrap();
    SdfReader::new(BufReader::new(f))
        .filter_map(|rec| rec.ok())
        .map(|rec| {
            let m = &rec.mol;
            let rb = ring_bundle(m);
            let mw = molecular_weight(m);
            let (logp, _) = logp_and_mr(m);
            DescriptorOutput { mw, logp, tpsa: tpsa(m), hbd: hbd_count(m), hba: rb.hba_count,
                lipinski_passes: mw <= 500.0 }
        })
        .collect()
}
```

## Why Rust for a chemistry API?

| Concern | Rust + chematic |
|---------|----------------|
| Throughput | ~3.6 µs/mol ECFP4 (5–14× faster than RDKit) |
| Memory safety | Enforced at compile time, no segfaults |
| Deployment | Single binary, no Python runtime |
| Embed in LIMS/ELN | `cargo add chematic`, no C++ toolchain |
| WASM companion | Same Rust code compiles to WASM for the browser |

## Related APIs

- `chematic_smiles::parse_smiles(&str)` — `→ Molecule`
- `chematic_chem::{molecular_weight, logp_and_mr, tpsa, hbd_count, ring_bundle}` — descriptor functions
- `chematic_fp::ecfp4(&mol)` / `tanimoto_bitvec(&fp1, &fp2)` — fingerprints + similarity
- `chematic_mol::SdfReader` — streaming SDF parser
- `chematic = { version = "1.0.4", features = ["chem", "fp", "smarts"] }` — minimal feature set for a descriptor API
