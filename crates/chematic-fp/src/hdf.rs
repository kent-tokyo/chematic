//! Hyper-Dimensional Fingerprints (HDF) for molecules.
//!
//! HDF encodes a molecule as a unit-norm float32 vector in a high-dimensional
//! space using Hyperdimensional Computing (HDC) operations.  Unlike hash-based
//! fingerprints (ECFP), there is no collision: every distinct atom environment
//! maps to a unique HD vector whose similarity can be measured with cosine dot.
//!
//! ## Algorithm
//!
//! 1. **Basis**: each (element, charge) pair → deterministic d-dim bipolar vector (±1).
//! 2. **Neighborhood**: for atom *i* at hop *h*, bind (element-wise multiply)
//!    the cyclic-shifted basis of *i* with the bases of *h*-hop neighbors.
//! 3. **Superposition**: sum all contributions across all atoms and hops.
//! 4. **Normalize** to unit sphere → `Vec<f32>` of length `dim`.
//!
//! Cosine similarity between two normalized HDF vectors equals their dot product,
//! making nearest-neighbor search trivial.
//!
//! ## References
//! - "Hyper-Dimensional Fingerprints as Molecular Representations" (arXiv 2026)
//! - Kanerva, P. "Hyperdimensional Computing" (Cognitive Computation 2009)

#![forbid(unsafe_code)]

use chematic_core::{AtomIdx, Molecule};

// ─── Config & Output ─────────────────────────────────────────────────────────

/// Configuration for HDF fingerprint generation.
#[derive(Debug, Clone)]
pub struct HdfConfig {
    /// Vector dimension. Higher → less collision, more memory. Default: 1024.
    pub dim: usize,
    /// Neighborhood radius (hops). Default: 2.
    pub radius: usize,
    /// Reproducibility seed. Default: 42.
    pub seed: u64,
}

impl Default for HdfConfig {
    fn default() -> Self {
        HdfConfig {
            dim: 1024,
            radius: 2,
            seed: 42,
        }
    }
}

/// A normalized HDF fingerprint vector (unit-norm f32 slice).
///
/// Use [`cosine_hdf`] for similarity between two HDF fingerprints of the same config.
#[derive(Debug, Clone)]
pub struct HdfFp(pub Vec<f32>);

// ─── Public API ──────────────────────────────────────────────────────────────

/// Compute an HDF fingerprint with default config (dim=1024, radius=2, seed=42).
pub fn hdf_default(mol: &Molecule) -> HdfFp {
    hdf(mol, &HdfConfig::default())
}

/// Compute an HDF fingerprint for `mol` using the given `config`.
pub fn hdf(mol: &Molecule, config: &HdfConfig) -> HdfFp {
    let d = config.dim;
    let n = mol.atom_count();

    if n == 0 {
        return HdfFp(vec![0.0f32; d]);
    }

    // Pre-compute per-atom basis vectors (±1 floats, deterministic)
    let bases: Vec<Vec<f32>> = (0..n)
        .map(|i| {
            let atom = mol.atom(AtomIdx(i as u32));
            random_bipolar(
                d,
                atom_seed(atom.element.atomic_number(), atom.charge, config.seed),
            )
        })
        .collect();

    let mut mol_vec = vec![0.0f32; d];

    for i in 0..n {
        // Hop 0: atom's own basis
        add_scaled(&mut mol_vec, &bases[i], 1.0);

        if config.radius == 0 {
            continue;
        }

        // BFS expansion up to radius
        let mut visited = vec![false; n];
        visited[i] = true;

        // frontier = immediate neighbors (hop 1)
        let mut frontier: Vec<usize> = mol
            .neighbors(AtomIdx(i as u32))
            .map(|(nb, _)| {
                let j = nb.0 as usize;
                visited[j] = true;
                j
            })
            .collect();

        // accumulated_binding starts as the center atom's basis
        let mut binding = bases[i].clone();

        for hop in 1..=config.radius {
            if frontier.is_empty() {
                break;
            }

            // Permute (cyclic-shift) the accumulated binding vector by `hop` positions
            // so each hop level produces a distinct encoding.
            let shifted = cyclic_shift(&binding, hop * 7 + 3); // prime-ish stride

            for &j in &frontier {
                // Bind: element-wise multiply shifted center with neighbor basis
                let bound = elementwise_product(&shifted, &bases[j]);
                add_scaled(&mut mol_vec, &bound, 1.0);
            }

            // Expand frontier one more hop
            let mut next: Vec<usize> = Vec::new();
            for &j in &frontier {
                for (nb, _) in mol.neighbors(AtomIdx(j as u32)) {
                    let k = nb.0 as usize;
                    if !visited[k] {
                        visited[k] = true;
                        next.push(k);
                    }
                }
            }

            // Accumulate the frontier atoms into binding for the next hop
            for &j in &frontier {
                elementwise_multiply_into(&mut binding, &bases[j]);
            }

            frontier = next;
        }
    }

    normalize(&mut mol_vec);
    HdfFp(mol_vec)
}

/// Cosine similarity between two HDF fingerprints (both must be normalized).
///
/// Returns a value in [-1, 1] where 1 = identical, 0 = orthogonal.
/// Panics in debug if the two vectors have different lengths.
pub fn cosine_hdf(a: &HdfFp, b: &HdfFp) -> f32 {
    debug_assert_eq!(a.0.len(), b.0.len(), "HdfFp dimensions must match");
    a.0.iter().zip(&b.0).map(|(x, y)| x * y).sum()
}

// ─── HD operations ───────────────────────────────────────────────────────────

/// Generate a d-dimensional bipolar (±1) vector seeded deterministically.
fn random_bipolar(d: usize, seed: u64) -> Vec<f32> {
    let mut state = if seed == 0 { 1 } else { seed };
    let mut out = Vec::with_capacity(d);
    while out.len() < d {
        // xorshift64
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        // Each u64 gives 64 bits
        let bits = state;
        for b in 0..64 {
            if out.len() == d {
                break;
            }
            out.push(if (bits >> b) & 1 == 0 {
                -1.0f32
            } else {
                1.0f32
            });
        }
    }
    out
}

/// Cyclic-shift a vector left by `k % d` positions.
fn cyclic_shift(v: &[f32], k: usize) -> Vec<f32> {
    let d = v.len();
    if d == 0 {
        return Vec::new();
    }
    let k = k % d;
    let mut out = Vec::with_capacity(d);
    out.extend_from_slice(&v[k..]);
    out.extend_from_slice(&v[..k]);
    out
}

/// Element-wise product (binding) of two equal-length vectors.
fn elementwise_product(a: &[f32], b: &[f32]) -> Vec<f32> {
    a.iter().zip(b).map(|(x, y)| x * y).collect()
}

/// Bind `binding` in-place: element-wise multiply by `factor`.
fn elementwise_multiply_into(binding: &mut [f32], factor: &[f32]) {
    for (b, &f) in binding.iter_mut().zip(factor) {
        *b *= f;
    }
}

/// Add `scale * src` into `dst`.
fn add_scaled(dst: &mut [f32], src: &[f32], scale: f32) {
    for (d, &s) in dst.iter_mut().zip(src) {
        *d += scale * s;
    }
}

/// Normalize `v` to unit length in-place; no-op if length is ~0.
fn normalize(v: &mut [f32]) {
    let norm: f32 = v.iter().map(|x| x * x).sum::<f32>().sqrt();
    if norm > 1e-9 {
        for x in v.iter_mut() {
            *x /= norm;
        }
    }
}

/// Deterministic seed for an atom given its element, charge, and global seed.
fn atom_seed(atomic_num: u8, charge: i8, global_seed: u64) -> u64 {
    global_seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(atomic_num as u64)
        .wrapping_add((charge as i64 + 128) as u64)
        .wrapping_mul(2654435761)
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn mol(smi: &str) -> chematic_core::Molecule {
        parse(smi).expect("parse SMILES")
    }

    fn hdf_default_mol(smi: &str) -> HdfFp {
        hdf_default(&mol(smi))
    }

    #[test]
    fn same_molecule_cosine_one() {
        let fp = hdf_default_mol("CCO");
        let fp2 = hdf_default_mol("CCO");
        let cos = cosine_hdf(&fp, &fp2);
        assert!((cos - 1.0).abs() < 1e-5, "same molecule cosine={cos}");
    }

    #[test]
    fn unit_norm() {
        for smi in ["CCO", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"] {
            let fp = hdf_default_mol(smi);
            let norm: f32 = fp.0.iter().map(|x| x * x).sum::<f32>().sqrt();
            assert!((norm - 1.0).abs() < 1e-5, "{smi}: norm={norm}");
        }
    }

    #[test]
    fn ethanol_vs_benzene_less_similar_than_propanol() {
        let eth = hdf_default_mol("CCO");
        let ben = hdf_default_mol("c1ccccc1");
        let pro = hdf_default_mol("CCCO");
        let cos_eb = cosine_hdf(&eth, &ben);
        let cos_ep = cosine_hdf(&eth, &pro);
        assert!(
            cos_ep > cos_eb,
            "propanol ({cos_ep:.3}) should be more similar to ethanol than benzene ({cos_eb:.3})"
        );
    }

    #[test]
    fn different_dims_work() {
        for dim in [64, 256, 512, 1024, 2048] {
            let config = HdfConfig {
                dim,
                radius: 2,
                seed: 42,
            };
            let fp = hdf(&mol("CCO"), &config);
            assert_eq!(fp.0.len(), dim);
            let norm: f32 = fp.0.iter().map(|x| x * x).sum::<f32>().sqrt();
            assert!((norm - 1.0).abs() < 1e-4, "dim={dim} norm={norm}");
        }
    }

    #[test]
    fn different_radius_gives_different_vectors() {
        let mol = mol("c1ccccc1");
        let r0 = hdf(
            &mol,
            &HdfConfig {
                dim: 256,
                radius: 0,
                seed: 1,
            },
        );
        let r2 = hdf(
            &mol,
            &HdfConfig {
                dim: 256,
                radius: 2,
                seed: 1,
            },
        );
        let cos = cosine_hdf(&r0, &r2);
        assert!(
            cos < 0.999,
            "radius-0 and radius-2 should differ; cosine={cos}"
        );
    }

    #[test]
    fn different_seeds_give_different_vectors() {
        let m = mol("CCO");
        let fp1 = hdf(
            &m,
            &HdfConfig {
                dim: 256,
                radius: 2,
                seed: 1,
            },
        );
        let fp2 = hdf(
            &m,
            &HdfConfig {
                dim: 256,
                radius: 2,
                seed: 2,
            },
        );
        let cos = cosine_hdf(&fp1, &fp2);
        assert!(
            cos < 0.999,
            "different seeds should give different vectors; cosine={cos}"
        );
    }

    #[test]
    fn charged_atom_differs_from_neutral() {
        // [NH4+] vs NH3: different charge → different basis → different HDF
        let charged = hdf_default_mol("[NH4+]");
        let neutral = hdf_default_mol("N");
        let cos = cosine_hdf(&charged, &neutral);
        assert!(cos < 0.99, "charged vs neutral too similar: cosine={cos}");
    }

    #[test]
    fn empty_molecule_returns_zero_vector() {
        let empty = chematic_core::MoleculeBuilder::new().build();
        let fp = hdf_default(&empty);
        assert!(fp.0.iter().all(|&x| x == 0.0));
    }
}
