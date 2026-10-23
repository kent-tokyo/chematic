# Embed bounded chemistry in a Rust service

## Problem

You are building a Rust service and need descriptors, similarity, or
substructure operations without adding a Python process or a C++ runtime.

## Solution

The `chematic` umbrella crate exposes selected chemistry modules behind feature
flags. The same process can validate input, run a supported operation, and
return a typed error. This guide is a small Axum sketch, not a production API
contract.

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
chematic = { version = "1.2.0", features = ["chem", "fp", "smarts"] }
axum      = "0.7"
serde     = { version = "1", features = ["derive"] }
tokio     = { version = "1", features = ["full"] }
```

## Descriptor endpoint

```rust
use axum::{extract::Json, http::StatusCode, routing::post, Router};
use chematic::{chem, smiles};
use serde::{Deserialize, Serialize};

#[derive(Deserialize)]
struct SmilesInput { smiles: String }

#[derive(Serialize)]
struct DescriptorOutput {
    mw: f64, logp: f64, tpsa: f64, hbd: usize, hba: usize,
    lipinski_passes: bool,
}

async fn descriptors(
    Json(body): Json<SmilesInput>,
) -> Result<Json<DescriptorOutput>, (StatusCode, String)> {
    let mol = smiles::parse(&body.smiles)
        .map_err(|error| (StatusCode::BAD_REQUEST, error.to_string()))?;
    let rb = chem::ring_bundle(&mol);
    let mw = chem::molecular_weight(&mol);
    let (logp, _) = chem::logp_and_mr(&mol);
    let tpsa = chem::tpsa(&mol);
    let hbd = chem::hbd_count(&mol);
    Ok(Json(DescriptorOutput {
        mw, logp, tpsa, hbd, hba: rb.hba_count,
        lipinski_passes: mw <= 500.0 && hbd <= 5 && rb.hba_count <= 10 && logp <= 5.0,
    }))
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
use axum::{extract::Json, http::StatusCode};
use chematic::{fp, smiles};

async fn similarity(
    Json(body): Json<SimilarityInput>,
) -> Result<Json<Vec<Hit>>, (StatusCode, String)> {
    let query_mol = smiles::parse(&body.query)
        .map_err(|error| (StatusCode::BAD_REQUEST, error.to_string()))?;
    let query_fp  = fp::ecfp4(&query_mol);

    let hits: Vec<Hit> = body.library
        .iter()
        .enumerate()
        .filter_map(|(i, smi)| {
            let mol = smiles::parse(smi).ok()?;
            let fp  = fp::ecfp4(&mol);
            let sim = fp::tanimoto_bitvec(&query_fp, &fp);
            (sim >= body.threshold).then(|| Hit { index: i, smiles: smi.clone(), score: sim })
        })
        .collect();

    Ok(Json(hits))
}
```

This compact similarity sketch omits malformed library entries. A production
endpoint should retain their original indices and return a separate error list
instead of silently changing the denominator.

## SDF batch processing

```rust
use chematic::{chem, mol::SdfFileReader};
use std::io::BufReader;

fn process_sdf(path: &str) -> Vec<DescriptorOutput> {
    let f = std::fs::File::open(path).unwrap();
    SdfFileReader::new(BufReader::new(f))
        .filter_map(|rec| rec.ok())
        .map(|rec| {
            let m = &rec.mol;
            let rb = chem::ring_bundle(m);
            let mw = chem::molecular_weight(m);
            let (logp, _) = chem::logp_and_mr(m);
            DescriptorOutput { mw, logp, tpsa: chem::tpsa(m), hbd: chem::hbd_count(m), hba: rb.hba_count,
                lipinski_passes: mw <= 500.0 }
        })
        .collect()
}
```

This compact SDF sketch also drops malformed records. In a production import,
return a row outcome containing the source index, status, and error so the
reported denominator still equals the number of submitted records.

## Why embed chematic here?

| Requirement | What chematic provides |
|---------|----------------|
| Deployment | Rust crates without a Python runtime or shared C++ chemistry library |
| Failure handling | `Result`-based parse and checked-operation errors |
| Resource policy | Documented limits for bounded operations and file readers |
| Browser companion | Related APIs are exposed from the same Rust core through WebAssembly |

Performance depends on the exact operation and build. Use the
[benchmark method](../reference/benchmark.md) and remeasure the endpoint, including
parsing, serialization, rejected rows, and peak memory.

## Related APIs

- `chematic_smiles::parse_smiles(&str)` — `→ Molecule`
- `chematic_chem::{molecular_weight, logp_and_mr, tpsa, hbd_count, ring_bundle}` — descriptor functions
- `chematic_fp::ecfp4(&mol)` / `tanimoto_bitvec(&fp1, &fp2)` — fingerprints + similarity
- `chematic_mol::SdfReader` — streaming SDF parser
- `chematic = { version = "1.2.0", features = ["chem", "fp", "smarts"] }` — minimal feature set for a descriptor API
