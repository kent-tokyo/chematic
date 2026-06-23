//! Spectrophores — 3D molecular fingerprints.
//!
//! Computes a 48-element vector (4 atomic properties × 12 probe positions)
//! that encodes the 3D electrostatic, lipophilic, aromatic, and H-bond
//! character of a molecule's surface.
//!
//! Originally developed by Silicos-it; patent expired 2024.
//!
//! # Algorithm
//!
//! 1. Compute per-heavy-atom values for four properties:
//!    - Electrostatic: Gasteiger partial charges
//!    - Lipophilic: Crippen LogP contributions
//!    - Aromatic: 1.0 for aromatic atoms, 0.0 otherwise
//!    - H-bond: +1.0 donor, −1.0 acceptor, 0.0 otherwise
//!
//! 2. Place 12 probe atoms at the vertices of a regular icosahedron
//!    scaled to `cage_radius = max_dist_from_centroid + cage_buffer`.
//!
//! 3. For each (probe, property) pair compute:
//!    `S = Σ_i  prop[i] / dist(probe, atom_i)^exponent`
//!
//! 4. Return the 48 values as `Vec<f64>` in row-major order:
//!    `[elec×12, lipo×12, arom×12, hbond×12]`.
//!
//! # Reference
//! Silicos-it Spectrophores, J. Cheminformatics (2018) 10:48.

use chematic_core::{AtomIdx, Molecule};
use chematic_perception::{FeatureType, detect_features};

use crate::coords::Coords3D;

// ---------------------------------------------------------------------------
// Regular icosahedron vertices (normalized to unit sphere)
// φ = (1 + √5) / 2,  |v| = √(1 + φ²) = √(2 + φ)
// ---------------------------------------------------------------------------

const PHI: f64 = 1.618_033_988_749_895; // (1 + √5) / 2
const NORM: f64 = 1.902_113_032_590_307; // √(2 + φ)

/// 12 vertices of a regular icosahedron, normalized to the unit sphere.
const ICOSAHEDRON_VERTICES: [[f64; 3]; 12] = [
    [0.0, 1.0 / NORM, PHI / NORM],
    [0.0, 1.0 / NORM, -PHI / NORM],
    [0.0, -1.0 / NORM, PHI / NORM],
    [0.0, -1.0 / NORM, -PHI / NORM],
    [1.0 / NORM, PHI / NORM, 0.0],
    [1.0 / NORM, -PHI / NORM, 0.0],
    [-1.0 / NORM, PHI / NORM, 0.0],
    [-1.0 / NORM, -PHI / NORM, 0.0],
    [PHI / NORM, 0.0, 1.0 / NORM],
    [PHI / NORM, 0.0, -1.0 / NORM],
    [-PHI / NORM, 0.0, 1.0 / NORM],
    [-PHI / NORM, 0.0, -1.0 / NORM],
];

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------

/// Normalization applied to the raw Spectrophores vector.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum SpectrophoresNorm {
    /// No normalization — return raw interaction energies.
    #[default]
    None,
    /// Per-property Z-score: subtract mean, divide by std dev (for each of 4 channels of 12).
    ZScore,
    /// Normalize the full 48-element vector to unit L2 norm.
    L2,
}

/// Configuration for Spectrophores computation.
#[derive(Debug, Clone)]
pub struct SpectrophoresConfig {
    /// Normalization applied to the output vector.
    pub normalize: SpectrophoresNorm,
    /// Extra buffer (Å) added to the maximum atom–centroid distance to set cage radius.
    pub cage_buffer: f64,
    /// Exponent for distance weighting: `1 / dist^exponent`.
    pub exponent: u32,
}

impl Default for SpectrophoresConfig {
    fn default() -> Self {
        Self {
            normalize: SpectrophoresNorm::None,
            cage_buffer: 3.0,
            exponent: 2,
        }
    }
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Compute the Spectrophores 3D fingerprint for a molecule.
///
/// Returns a 48-element `Vec<f64>` — four blocks of 12 probe values for the
/// electrostatic, lipophilic, aromatic, and H-bond properties respectively.
///
/// `coords` must supply one `Point3` per heavy atom in the same atom order as
/// `mol`.  Use [`crate::generate_coords`] or [`crate::generate_and_minimize_constrained`]
/// to obtain 3D coordinates.
///
/// # Example (Rust)
/// ```ignore
/// let mol = chematic_smiles::parse("c1ccccc1").unwrap();
/// let coords = chematic_3d::generate_and_minimize_constrained(&mol);
/// let fp = chematic_3d::spectrophores(&mol, &coords, &Default::default());
/// assert_eq!(fp.len(), 48);
/// ```
pub fn spectrophores(mol: &Molecule, coords: &Coords3D, config: &SpectrophoresConfig) -> Vec<f64> {
    // Collect heavy atoms and their positions.
    let heavy: Vec<AtomIdx> = mol
        .atoms()
        .filter(|(_, a)| a.element.atomic_number() != 1)
        .map(|(idx, _)| idx)
        .collect();

    if heavy.is_empty() {
        return vec![0.0; 48];
    }

    let positions: Vec<[f64; 3]> = heavy
        .iter()
        .map(|&idx| {
            let p = coords.get(idx);
            [p.x, p.y, p.z]
        })
        .collect();

    // Compute 4 property arrays.
    let elec = electrostatic_props(mol, &heavy);
    let lipo = lipophilic_props(mol, &heavy);
    let arom = aromatic_props(mol, &heavy);
    let hbond = hbond_props(mol, &heavy);
    let props: [&[f64]; 4] = [&elec, &lipo, &arom, &hbond];

    // Translate centroid to origin.
    let cx = positions.iter().map(|p| p[0]).sum::<f64>() / positions.len() as f64;
    let cy = positions.iter().map(|p| p[1]).sum::<f64>() / positions.len() as f64;
    let cz = positions.iter().map(|p| p[2]).sum::<f64>() / positions.len() as f64;
    let centered: Vec<[f64; 3]> = positions
        .iter()
        .map(|p| [p[0] - cx, p[1] - cy, p[2] - cz])
        .collect();

    // Cage radius: max atom distance from centroid + buffer.
    let max_dist = centered
        .iter()
        .map(|p| (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt())
        .fold(0.0_f64, f64::max);
    let cage_r = max_dist + config.cage_buffer;

    // Build 48 values: for each property × probe.
    let mut result = Vec::with_capacity(48);
    for prop_vals in props {
        for vertex in &ICOSAHEDRON_VERTICES {
            let probe = [vertex[0] * cage_r, vertex[1] * cage_r, vertex[2] * cage_r];
            let s: f64 = centered
                .iter()
                .zip(prop_vals.iter())
                .map(|(atom_pos, &q)| {
                    let dx = probe[0] - atom_pos[0];
                    let dy = probe[1] - atom_pos[1];
                    let dz = probe[2] - atom_pos[2];
                    let d2 = dx * dx + dy * dy + dz * dz;
                    if d2 < 1e-12 {
                        return 0.0; // guard against degenerate case
                    }
                    let d = d2.sqrt();
                    q / d.powi(config.exponent as i32)
                })
                .sum();
            result.push(s);
        }
    }

    apply_normalization(&mut result, config.normalize);
    result
}

/// Tanimoto-like similarity between two Spectrophores vectors.
///
/// Uses the USR formula: `S = 1 / (1 + mean|a − b|)` giving values in (0, 1].
/// Returns 1.0 for identical vectors, approaches 0 for very different molecules.
pub fn tanimoto_spectrophores(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(
        a.len(),
        b.len(),
        "Spectrophores vectors must have equal length"
    );
    let n = a.len() as f64;
    if n == 0.0 {
        return 1.0;
    }
    let mean_abs_diff = a
        .iter()
        .zip(b.iter())
        .map(|(x, y)| (x - y).abs())
        .sum::<f64>()
        / n;
    1.0 / (1.0 + mean_abs_diff)
}

// ---------------------------------------------------------------------------
// Property computation
// ---------------------------------------------------------------------------

fn electrostatic_props(mol: &Molecule, heavy: &[AtomIdx]) -> Vec<f64> {
    let all_charges = chematic_chem::gasteiger_charges(mol);
    heavy
        .iter()
        .map(|&idx| *all_charges.get(idx.0 as usize).unwrap_or(&0.0))
        .collect()
}

fn lipophilic_props(mol: &Molecule, heavy: &[AtomIdx]) -> Vec<f64> {
    let all_logp = chematic_chem::logp_crippen_per_atom(mol);
    heavy
        .iter()
        .map(|&idx| *all_logp.get(idx.0 as usize).unwrap_or(&0.0))
        .collect()
}

fn aromatic_props(mol: &Molecule, heavy: &[AtomIdx]) -> Vec<f64> {
    heavy
        .iter()
        .map(|&idx| if mol.atom(idx).aromatic { 1.0 } else { 0.0 })
        .collect()
}

fn hbond_props(mol: &Molecule, heavy: &[AtomIdx]) -> Vec<f64> {
    let n = mol.atom_count();
    let mut scores = vec![0.0f64; n];

    for feat in detect_features(mol) {
        let contribution = match feat.ftype {
            FeatureType::Donor => 1.0,
            FeatureType::Acceptor => -1.0,
            _ => continue,
        };
        let i = feat.atom.0 as usize;
        if i < n {
            scores[i] += contribution;
        }
    }

    heavy.iter().map(|&idx| scores[idx.0 as usize]).collect()
}

// ---------------------------------------------------------------------------
// Normalization
// ---------------------------------------------------------------------------

fn apply_normalization(v: &mut [f64], mode: SpectrophoresNorm) {
    match mode {
        SpectrophoresNorm::None => {}
        SpectrophoresNorm::ZScore => {
            // Normalize each of 4 property channels (12 values each) independently.
            for chunk in v.chunks_mut(12) {
                let n = chunk.len() as f64;
                let mean = chunk.iter().sum::<f64>() / n;
                let var = chunk.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n;
                let std = var.sqrt();
                if std > 1e-12 {
                    for x in chunk.iter_mut() {
                        *x = (*x - mean) / std;
                    }
                } else {
                    for x in chunk.iter_mut() {
                        *x -= mean;
                    }
                }
            }
        }
        SpectrophoresNorm::L2 => {
            let norm = v.iter().map(|x| x * x).sum::<f64>().sqrt();
            if norm > 1e-12 {
                for x in v.iter_mut() {
                    *x /= norm;
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn make_coords(mol: &Molecule) -> Coords3D {
        crate::dg::generate_coords(mol)
    }

    #[test]
    fn benzene_48_elements_finite() {
        let mol = parse("c1ccccc1").unwrap();
        let coords = make_coords(&mol);
        let fp = spectrophores(&mol, &coords, &Default::default());
        assert_eq!(fp.len(), 48, "must return 48 elements");
        assert!(
            fp.iter().all(|v| v.is_finite()),
            "all values must be finite"
        );
    }

    #[test]
    fn empty_molecule_returns_zeros() {
        let mol = parse("[H][H]").unwrap();
        let coords = make_coords(&mol);
        let fp = spectrophores(&mol, &coords, &Default::default());
        assert_eq!(fp.len(), 48);
        assert!(fp.iter().all(|v| *v == 0.0));
    }

    #[test]
    fn same_molecule_tanimoto_is_one() {
        let mol = parse("CC(=O)O").unwrap();
        let coords = make_coords(&mol);
        let fp = spectrophores(&mol, &coords, &Default::default());
        let sim = tanimoto_spectrophores(&fp, &fp);
        assert!((sim - 1.0).abs() < 1e-10, "tanimoto of self must be 1.0");
    }

    #[test]
    fn zscore_normalization_finite() {
        let mol = parse("CN1C=NC2=C1C(=O)N(C)C(=O)N2C").unwrap(); // caffeine
        let coords = make_coords(&mol);
        let config = SpectrophoresConfig {
            normalize: SpectrophoresNorm::ZScore,
            ..Default::default()
        };
        let fp = spectrophores(&mol, &coords, &config);
        assert_eq!(fp.len(), 48);
        assert!(fp.iter().all(|v| v.is_finite()));
    }

    #[test]
    fn different_molecules_tanimoto_less_than_one() {
        let benzene = parse("c1ccccc1").unwrap();
        let water = parse("O").unwrap();
        let c1 = make_coords(&benzene);
        let c2 = make_coords(&water);
        let fp1 = spectrophores(&benzene, &c1, &Default::default());
        let fp2 = spectrophores(&water, &c2, &Default::default());
        let sim = tanimoto_spectrophores(&fp1, &fp2);
        assert!(sim < 1.0, "benzene and water should differ: {sim}");
    }
}
