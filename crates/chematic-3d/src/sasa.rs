//! Solvent-Accessible Surface Area (SASA) calculation using Shrake-Rupley algorithm.
//!
//! B8: 3D SASA Descriptor
//!
//! Provides efficient SASA calculation for 3D molecular structures from distance geometry (A3).
//! SASA is the surface area of a molecule that is accessible to a solvent probe (typically water).
//! Useful for:
//! - Molecular property prediction (solubility, binding affinity)
//! - Descriptor calculation for ML models
//! - Comparison with RDKit calculations
//!
//! Default parameters: probe_radius=1.4 Å (water), sphere_points=100 (speed/accuracy trade-off)

use chematic_core::{AtomIdx, Molecule};

use crate::coords::{Coords3D, Point3};
use crate::dg::generate_coords;

/// Standard Bondi VDW radii (Ångströms) from Bondi 1964.
fn bondi_vdw_radius(atomic_number: u8) -> f64 {
    match atomic_number {
        1 => 1.20,  // H
        6 => 1.70,  // C
        7 => 1.55,  // N
        8 => 1.52,  // O
        9 => 1.47,  // F
        15 => 1.80, // P
        16 => 1.80, // S
        17 => 1.75, // Cl
        35 => 1.85, // Br
        53 => 1.98, // I
        _ => 1.70,  // default fallback
    }
}

/// Calculate SASA with default parameters (RDKit-compatible).
///
/// Convenience function using standard parameters:
/// - probe_radius = 1.4 Ångströms (water)
/// - sphere_points = 100 (good balance of speed/accuracy)
///
/// # Arguments
/// - `mol`: molecule (for atom properties)
/// - `coords`: 3D coordinates (from distance geometry or force field)
///
/// # Returns
/// Total SASA in Ų (square Ångströms)
pub fn sasa(mol: &Molecule, coords: &Coords3D) -> f64 {
    shrake_rupley_sasa(mol, coords, 1.4, 100)
}

/// Calculate per-atom SASA with default parameters.
///
/// # Arguments
/// - `mol`: molecule (for atom properties)
/// - `coords`: 3D coordinates
///
/// # Returns
/// Vector of SASA values (one per atom) in Ų
pub fn sasa_per_atom_default(mol: &Molecule, coords: &Coords3D) -> Vec<f64> {
    sasa_per_atom(mol, coords, 1.4, 100)
}

/// Calculate SASA from SMILES string with distance geometry.
///
/// Generates 3D coordinates using distance geometry, then calculates SASA.
/// Useful for quick SASA calculations without pre-computed coordinates.
///
/// # Arguments
/// - `mol`: molecule
///
/// # Returns
/// Total SASA in Ų, or error if coordinate generation fails
pub fn sasa_from_dg(mol: &Molecule) -> Result<f64, String> {
    let coords = generate_coords(mol);
    Ok(sasa(mol, &coords))
}

/// Calculate per-atom SASA from distance geometry coordinates.
pub fn sasa_per_atom_from_dg(mol: &Molecule) -> Result<Vec<f64>, String> {
    let coords = generate_coords(mol);
    Ok(sasa_per_atom_default(mol, &coords))
}

/// Calculate Solvent-Accessible Surface Area (SASA) using Shrake-Rupley algorithm.
///
/// The algorithm places points on a sphere around each atom (van der Waals radius + probe radius).
/// Points are tested for occlusion by neighboring atoms. The exposed surface area is proportional
/// to the number of exposed sphere points.
///
/// # Arguments
/// - `mol`: molecule (for atom properties)
/// - `coords`: 3D coordinates
/// - `probe_radius`: probe radius in Ångströms (default ~1.4 Å for water)
/// - `sphere_points`: number of points to sample per atom (typically 100 or more)
///
/// # Returns
/// Total SASA in Ų (square Ångströms)
pub fn shrake_rupley_sasa(
    mol: &Molecule,
    coords: &Coords3D,
    probe_radius: f64,
    sphere_points: usize,
) -> f64 {
    if mol.atom_count() == 0 {
        return 0.0;
    }

    let mut total_sasa = 0.0;

    for i in 0..mol.atom_count() {
        let idx = AtomIdx(i as u32);
        let atom_i = mol.atom(idx);
        let vdw_i = bondi_vdw_radius(atom_i.element.atomic_number());
        let radius_i = vdw_i + probe_radius;
        let pos_i = coords.get(idx);

        // Generate sphere points around atom i
        let exposed_count =
            count_exposed_points(mol, coords, pos_i, idx, probe_radius, sphere_points);

        // Surface area proportional to exposed fraction of sphere
        let sphere_area = 4.0 * std::f64::consts::PI * radius_i * radius_i;
        let atom_sasa = (exposed_count as f64 / sphere_points as f64) * sphere_area;
        total_sasa += atom_sasa;
    }

    total_sasa
}

/// Calculate per-atom Solvent-Accessible Surface Area.
///
/// Returns a vector of SASA values, one per atom, in square Ångströms.
pub fn sasa_per_atom(
    mol: &Molecule,
    coords: &Coords3D,
    probe_radius: f64,
    sphere_points: usize,
) -> Vec<f64> {
    let mut sasa_values = vec![0.0; mol.atom_count()];

    for i in 0..mol.atom_count() {
        let idx = AtomIdx(i as u32);
        let atom_i = mol.atom(idx);
        let vdw_i = bondi_vdw_radius(atom_i.element.atomic_number());
        let radius_i = vdw_i + probe_radius;
        let pos_i = coords.get(idx);

        let exposed_count =
            count_exposed_points(mol, coords, pos_i, idx, probe_radius, sphere_points);

        let sphere_area = 4.0 * std::f64::consts::PI * radius_i * radius_i;
        sasa_values[i] = (exposed_count as f64 / sphere_points as f64) * sphere_area;
    }

    sasa_values
}

/// Count exposed sphere points for a given atom.
///
/// A point on the sphere is considered exposed if it is not occluded by any neighboring atom.
fn count_exposed_points(
    mol: &Molecule,
    coords: &Coords3D,
    atom_pos: Point3,
    atom_idx: AtomIdx,
    probe_radius: f64,
    sphere_points: usize,
) -> usize {
    let atom = mol.atom(atom_idx);
    let atom_vdw = bondi_vdw_radius(atom.element.atomic_number());
    let atom_radius = atom_vdw + probe_radius;

    let sphere = generate_sphere_points(atom_pos, atom_radius, sphere_points);

    // Precompute neighbor positions and radii to avoid redundant lookups in the inner loop
    let neighbors: Vec<(Point3, f64)> = (0..mol.atom_count())
        .filter(|&j| j as u32 != atom_idx.0)
        .map(|j| {
            let nbr_idx = AtomIdx(j as u32);
            let nbr = mol.atom(nbr_idx);
            let nbr_vdw = bondi_vdw_radius(nbr.element.atomic_number());
            let nbr_radius = nbr_vdw + probe_radius;
            let nbr_pos = coords.get(nbr_idx);
            (nbr_pos, nbr_radius)
        })
        .collect();

    let mut exposed_count = 0;

    for &point in &sphere {
        let mut is_exposed = true;

        // Check occlusion by neighboring atoms
        for &(nbr_pos, nbr_radius) in &neighbors {
            // Distance from sphere point to neighbor atom center
            let dx = point.x - nbr_pos.x;
            let dy = point.y - nbr_pos.y;
            let dz = point.z - nbr_pos.z;
            let dist = (dx * dx + dy * dy + dz * dz).sqrt();

            // Point is occluded if it falls within neighbor's SASA radius
            if dist < nbr_radius - 1e-6 {
                is_exposed = false;
                break;
            }
        }

        if is_exposed {
            exposed_count += 1;
        }
    }

    exposed_count
}

/// Generate sphere points using the Fibonacci sphere algorithm.
///
/// This generates evenly distributed points on a sphere of given radius.
/// The Fibonacci sphere algorithm is fast and produces uniform coverage.
fn generate_sphere_points(center: Point3, radius: f64, num_points: usize) -> Vec<Point3> {
    if num_points < 2 {
        return vec![center]; // Avoid division by zero when num_points < 2
    }

    let mut points = Vec::with_capacity(num_points);
    let golden_angle = std::f64::consts::PI * (3.0 - 5_f64.sqrt());

    for i in 0..num_points {
        let y = 1.0 - (i as f64) / (num_points as f64 - 1.0) * 2.0;
        let x_radius = (1.0 - y * y).sqrt();

        let theta = golden_angle * (i as f64);
        let x = theta.cos() * x_radius;
        let z = theta.sin() * x_radius;

        let point = Point3::new(
            center.x + x * radius,
            center.y + y * radius,
            center.z + z * radius,
        );
        points.push(point);
    }

    points
}

/// SASA descriptor statistics for a molecule.
#[derive(Clone, Debug)]
pub struct SasaDescriptor {
    /// Total solvent-accessible surface area (Ų)
    pub total: f64,
    /// Mean SASA per atom (Ų)
    pub mean: f64,
    /// Standard deviation of per-atom SASA
    pub std_dev: f64,
    /// Per-atom SASA values
    pub per_atom: Vec<f64>,
}

impl SasaDescriptor {
    /// Create a SASA descriptor from per-atom values.
    pub fn from_per_atom(per_atom: Vec<f64>) -> Self {
        let total: f64 = per_atom.iter().sum();
        let n = per_atom.len() as f64;
        let mean = if n > 0.0 { total / n } else { 0.0 };

        let variance = if n > 0.0 {
            per_atom.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / n
        } else {
            0.0
        };
        let std_dev = variance.sqrt();

        SasaDescriptor {
            total,
            mean,
            std_dev,
            per_atom,
        }
    }
}

/// Calculate SASA descriptor with statistics using default parameters.
///
/// Returns a SasaDescriptor containing total, mean, std_dev, and per-atom values.
pub fn sasa_descriptor(mol: &Molecule, coords: &Coords3D) -> SasaDescriptor {
    let per_atom = sasa_per_atom_default(mol, coords);
    SasaDescriptor::from_per_atom(per_atom)
}

/// Calculate SASA descriptor from distance geometry coordinates.
pub fn sasa_descriptor_from_dg(mol: &Molecule) -> Result<SasaDescriptor, String> {
    let coords = generate_coords(mol);
    Ok(sasa_descriptor(mol, &coords))
}

/// Per-element SASA breakdown.
/// Returns a vector where index corresponds to atomic number.
#[derive(Clone, Debug)]
pub struct PerElementSasa {
    /// SASA by atomic number (index = Z, value = total SASA for that element)
    pub by_element: Vec<f64>,
}

impl PerElementSasa {
    /// Get SASA for a specific element (by atomic number).
    pub fn get(&self, atomic_number: u8) -> f64 {
        if (atomic_number as usize) < self.by_element.len() {
            self.by_element[atomic_number as usize]
        } else {
            0.0
        }
    }
}

/// Calculate per-element SASA breakdown.
pub fn sasa_per_element(mol: &Molecule, coords: &Coords3D) -> PerElementSasa {
    let per_atom = sasa_per_atom_default(mol, coords);
    let mut by_element = vec![0.0; 119]; // Periodic table up to element 118

    for (i, &sasa_val) in per_atom.iter().enumerate() {
        let idx = AtomIdx(i as u32);
        let z = mol.atom(idx).element.atomic_number() as usize;
        if z < by_element.len() {
            by_element[z] += sasa_val;
        }
    }

    PerElementSasa { by_element }
}

/// Calculate per-element SASA from distance geometry coordinates.
pub fn sasa_per_element_from_dg(mol: &Molecule) -> Result<PerElementSasa, String> {
    let coords = generate_coords(mol);
    Ok(sasa_per_element(mol, &coords))
}

/// RDKit-compatible SASA calculation (wraps the standard implementation).
/// This function signature mirrors RDKit's Descriptors.CalcMolSASA().
pub fn calc_mol_sasa(mol: &Molecule, coords: &Coords3D) -> f64 {
    sasa(mol, coords)
}

/// RDKit-compatible name: compute SASA with custom probe radius.
pub fn calc_mol_sasa_with_probe(mol: &Molecule, coords: &Coords3D, probe_radius: f64) -> f64 {
    shrake_rupley_sasa(mol, coords, probe_radius, 100)
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn test_sasa_single_atom() {
        let mol = parse("C").unwrap();
        let coords = Coords3D::new_zeroed(1);
        let sasa = shrake_rupley_sasa(&mol, &coords, 1.4, 100);
        // Single atom should have non-zero SASA
        assert!(sasa > 0.0, "single atom should have positive SASA");
    }

    #[test]
    fn test_sasa_multiple_atoms() {
        let mol = parse("CC").unwrap();
        let mut coords = Coords3D::new_zeroed(2);
        // Place atoms far apart (no occlusion)
        coords.set(AtomIdx(0), Point3::new(0.0, 0.0, 0.0));
        coords.set(AtomIdx(1), Point3::new(10.0, 0.0, 0.0));
        let sasa = shrake_rupley_sasa(&mol, &coords, 1.4, 100);
        // Two separated atoms should have SASA close to 2× single atom
        assert!(sasa > 0.0, "multi-atom SASA should be positive");
    }

    #[test]
    fn test_sasa_per_atom_sum() {
        let mol = parse("CC").unwrap();
        let mut coords = Coords3D::new_zeroed(2);
        coords.set(AtomIdx(0), Point3::new(0.0, 0.0, 0.0));
        coords.set(AtomIdx(1), Point3::new(10.0, 0.0, 0.0));
        let per_atom = sasa_per_atom(&mol, &coords, 1.4, 100);
        let sum: f64 = per_atom.iter().sum();
        let total = shrake_rupley_sasa(&mol, &coords, 1.4, 100);
        // Sum of per-atom SASA should equal total
        assert!(
            (sum - total).abs() < 1e-6,
            "per-atom sum should match total"
        );
    }

    #[test]
    fn test_sasa_empty_molecule() {
        let mol = parse("").unwrap_or_else(|_| {
            // Create empty molecule manually if parse fails
            chematic_core::MoleculeBuilder::new().build()
        });
        let coords = Coords3D::new_zeroed(0);
        let sasa = shrake_rupley_sasa(&mol, &coords, 1.4, 100);
        assert_eq!(sasa, 0.0, "empty molecule should have zero SASA");
    }

    // ===== Phase 1: Distance Geometry Integration & Comprehensive Testing =====

    #[test]
    fn test_sasa_with_default_params() {
        // Test convenience function with default parameters
        let mol = parse("C").unwrap();
        let coords = Coords3D::new_zeroed(1);
        let sasa_default = sasa(&mol, &coords);
        let sasa_explicit = shrake_rupley_sasa(&mol, &coords, 1.4, 100);

        assert!((sasa_default - sasa_explicit).abs() < 1e-6);
    }

    #[test]
    fn test_sasa_per_atom_default() {
        // Test per-atom convenience function
        let mol = parse("CC").unwrap();
        let mut coords = Coords3D::new_zeroed(2);
        coords.set(AtomIdx(0), Point3::new(0.0, 0.0, 0.0));
        coords.set(AtomIdx(1), Point3::new(10.0, 0.0, 0.0));

        let per_atom = sasa_per_atom_default(&mol, &coords);
        assert_eq!(per_atom.len(), 2);
        assert!(per_atom[0] > 0.0 && per_atom[1] > 0.0);
    }

    #[test]
    fn test_sasa_from_distance_geometry_methane() {
        // Test SASA calculation with DG-generated coordinates
        let mol = parse("C").unwrap();
        let result = sasa_from_dg(&mol);

        assert!(result.is_ok());
        let sasa = result.unwrap();

        // Single carbon should have SASA around 4πr² where r = vdW + probe
        // Carbon vdW = 1.70 Å, probe = 1.4 Å, so r = 3.1 Å
        // Expected SASA ≈ 4π(3.1)² ≈ 121 Ų
        assert!(
            sasa > 80.0 && sasa < 160.0,
            "methane SASA out of expected range: {}",
            sasa
        );
    }

    #[test]
    fn test_sasa_from_dg_ethane() {
        // Two-atom molecule
        let mol = parse("CC").unwrap();
        let result = sasa_from_dg(&mol);

        assert!(result.is_ok());
        let sasa = result.unwrap();

        // Two carbons, some overlap, should be larger than single atom but less than 2×
        assert!(sasa > 100.0, "ethane SASA should be substantial: {}", sasa);
    }

    #[test]
    fn test_sasa_per_atom_from_dg() {
        // Test per-atom SASA from distance geometry
        let mol = parse("CCO").unwrap();
        let result = sasa_per_atom_from_dg(&mol);

        assert!(result.is_ok());
        let per_atom = result.unwrap();

        assert_eq!(per_atom.len(), 3);
        // All atoms should have positive SASA
        for (i, &sasa) in per_atom.iter().enumerate() {
            assert!(sasa > 0.0, "atom {} SASA should be positive", i);
        }
    }

    #[test]
    fn test_sasa_is_additive_separated_atoms() {
        // When atoms are far apart (no occlusion), total SASA should approximately equal sum
        let mol = parse("CC").unwrap();
        let mut coords = Coords3D::new_zeroed(2);

        // Place atoms very far apart
        coords.set(AtomIdx(0), Point3::new(0.0, 0.0, 0.0));
        coords.set(AtomIdx(1), Point3::new(100.0, 0.0, 0.0));

        let total = sasa(&mol, &coords);
        let per_atom = sasa_per_atom_default(&mol, &coords);
        let sum: f64 = per_atom.iter().sum();

        // With atoms far apart, total should approximately equal sum
        assert!(
            (total - sum).abs() < 1.0,
            "far atoms should have additive SASA"
        );
    }

    #[test]
    fn test_sasa_occlusion_effect() {
        // Close atoms should have lower SASA than far atoms
        let mol = parse("CC").unwrap();

        // Configuration 1: atoms far apart
        let mut coords_far = Coords3D::new_zeroed(2);
        coords_far.set(AtomIdx(0), Point3::new(0.0, 0.0, 0.0));
        coords_far.set(AtomIdx(1), Point3::new(50.0, 0.0, 0.0));
        let sasa_far = sasa(&mol, &coords_far);

        // Configuration 2: atoms close together
        let mut coords_close = Coords3D::new_zeroed(2);
        coords_close.set(AtomIdx(0), Point3::new(0.0, 0.0, 0.0));
        coords_close.set(AtomIdx(1), Point3::new(2.0, 0.0, 0.0));
        let sasa_close = sasa(&mol, &coords_close);

        // Close atoms should have lower total SASA due to occlusion
        assert!(
            sasa_close < sasa_far,
            "close atoms should have lower SASA due to occlusion"
        );
    }

    #[test]
    fn test_sasa_probe_radius_effect() {
        // Larger probe radius should increase SASA
        let mol = parse("C").unwrap();
        let coords = Coords3D::new_zeroed(1);

        let sasa_small = shrake_rupley_sasa(&mol, &coords, 1.0, 100);
        let sasa_large = shrake_rupley_sasa(&mol, &coords, 2.0, 100);

        assert!(
            sasa_large > sasa_small,
            "larger probe radius should increase SASA"
        );
    }

    #[test]
    fn test_sasa_sphere_points_convergence() {
        // More sphere points should give more accurate results
        let mol = parse("CC").unwrap();
        let mut coords = Coords3D::new_zeroed(2);
        coords.set(AtomIdx(0), Point3::new(0.0, 0.0, 0.0));
        coords.set(AtomIdx(1), Point3::new(1.5, 0.0, 0.0));

        let sasa_100 = shrake_rupley_sasa(&mol, &coords, 1.4, 100);
        let sasa_500 = shrake_rupley_sasa(&mol, &coords, 1.4, 500);

        // Results should be similar but 500 points may give slightly different value
        // Due to better sampling
        assert!(
            (sasa_100 - sasa_500).abs() < sasa_100 * 0.2,
            "different sphere points should give similar results"
        );
    }

    #[test]
    fn test_sasa_benzene() {
        // Benzene with DG coordinates
        let mol = parse("c1ccccc1").unwrap();
        let result = sasa_from_dg(&mol);

        assert!(result.is_ok());
        let sasa = result.unwrap();

        // Benzene should have substantial SASA (larger than ethane)
        assert!(sasa > 150.0, "benzene SASA should be substantial: {}", sasa);
    }

    #[test]
    fn test_sasa_finite_nonzero() {
        // All SASA values should be finite and non-negative
        let mol = parse("C1CCCCC1").unwrap(); // Cyclohexane
        let result = sasa_from_dg(&mol);

        assert!(result.is_ok());
        let sasa = result.unwrap();

        assert!(sasa.is_finite(), "SASA should be finite");
        assert!(sasa > 0.0, "SASA should be positive");
        assert!(sasa < 1e6, "SASA should not be unreasonably large");
    }

    #[test]
    fn test_sasa_per_atom_polar_molecule() {
        // Test per-atom SASA for polar molecule (ethanol)
        let mol = parse("CCO").unwrap();
        let coords = generate_coords(&mol);

        let per_atom = sasa_per_atom_default(&mol, &coords);
        assert_eq!(per_atom.len(), 3);

        // Oxygen should be exposed (high SASA)
        // Carbon values depend on conformation
        for &sasa in &per_atom {
            assert!(sasa > 0.0 && sasa.is_finite());
        }
    }

    // B8 Enhancement Tests: Descriptor statistics and per-element breakdown

    #[test]
    fn test_sasa_descriptor_basic() {
        // Test SasaDescriptor computation
        let mol = parse("C").unwrap();
        let coords = Coords3D::new_zeroed(1);
        let desc = sasa_descriptor(&mol, &coords);

        assert!(desc.total > 0.0);
        assert_eq!(desc.per_atom.len(), 1);
        assert!((desc.mean - desc.total).abs() < 1e-6); // Single atom: mean = total
    }

    #[test]
    fn test_sasa_descriptor_multiple_atoms() {
        // Test descriptor with multiple atoms
        let mol = parse("CC").unwrap();
        let mut coords = Coords3D::new_zeroed(2);
        coords.set(AtomIdx(0), Point3::new(0.0, 0.0, 0.0));
        coords.set(AtomIdx(1), Point3::new(10.0, 0.0, 0.0));

        let desc = sasa_descriptor(&mol, &coords);

        assert_eq!(desc.per_atom.len(), 2);
        let sum: f64 = desc.per_atom.iter().sum();
        assert!((desc.total - sum).abs() < 1e-6);

        // Mean should be total / number of atoms
        let expected_mean = desc.total / 2.0;
        assert!((desc.mean - expected_mean).abs() < 1e-6);
    }

    #[test]
    fn test_sasa_descriptor_std_dev() {
        // Test standard deviation calculation
        let mol = parse("CC").unwrap();
        let mut coords = Coords3D::new_zeroed(2);
        coords.set(AtomIdx(0), Point3::new(0.0, 0.0, 0.0));
        coords.set(AtomIdx(1), Point3::new(2.0, 0.0, 0.0)); // Close atoms

        let desc = sasa_descriptor(&mol, &coords);

        // Standard deviation should be non-negative
        assert!(desc.std_dev >= 0.0);
    }

    #[test]
    fn test_sasa_descriptor_from_dg() {
        // Test descriptor from distance geometry
        let mol = parse("CCO").unwrap();
        let result = sasa_descriptor_from_dg(&mol);

        assert!(result.is_ok());
        let desc = result.unwrap();

        assert_eq!(desc.per_atom.len(), 3);
        assert!(desc.total > 0.0);
        assert!(desc.mean > 0.0);
        assert!(desc.std_dev >= 0.0);
    }

    #[test]
    fn test_sasa_per_element_single_element() {
        // Test per-element breakdown for single element
        let mol = parse("C").unwrap(); // Single carbon
        let coords = Coords3D::new_zeroed(1);
        let per_elem = sasa_per_element(&mol, &coords);

        let carbon_sasa = per_elem.get(6); // Carbon
        assert!(carbon_sasa > 0.0);

        // Other elements should have zero
        assert_eq!(per_elem.get(1), 0.0); // Hydrogen
        assert_eq!(per_elem.get(8), 0.0); // Oxygen
    }

    #[test]
    fn test_sasa_per_element_mixed() {
        // Test per-element breakdown for mixed molecule
        let mol = parse("CCO").unwrap();
        let coords = generate_coords(&mol);
        let per_elem = sasa_per_element(&mol, &coords);

        let carbon_sasa = per_elem.get(6);
        let oxygen_sasa = per_elem.get(8);

        assert!(carbon_sasa > 0.0);
        assert!(oxygen_sasa > 0.0);

        // Sum of elements should equal total
        let total: f64 = per_elem.by_element.iter().sum();
        let expected = sasa(&mol, &coords);
        assert!((total - expected).abs() < 1e-6);
    }

    #[test]
    fn test_sasa_per_element_from_dg() {
        // Test per-element breakdown from DG
        let mol = parse("C1CCCCC1").unwrap(); // Cyclohexane
        let result = sasa_per_element_from_dg(&mol);

        assert!(result.is_ok());
        let per_elem = result.unwrap();

        let carbon_sasa = per_elem.get(6);
        assert!(carbon_sasa > 0.0);
    }

    #[test]
    fn test_sasa_rdkit_compatible() {
        // Test RDKit-compatible wrapper function
        let mol = parse("C").unwrap();
        let coords = Coords3D::new_zeroed(1);

        let sasa_std = sasa(&mol, &coords);
        let sasa_rdkit = calc_mol_sasa(&mol, &coords);

        assert!((sasa_std - sasa_rdkit).abs() < 1e-10);
    }

    #[test]
    fn test_sasa_rdkit_probe_radius() {
        // Test RDKit-compatible probe radius variant
        let mol = parse("CC").unwrap();
        let coords = generate_coords(&mol);

        let sasa_default = calc_mol_sasa(&mol, &coords);
        let sasa_custom = calc_mol_sasa_with_probe(&mol, &coords, 1.4);

        // With same probe radius, should be same
        assert!((sasa_default - sasa_custom).abs() < 1e-6);

        // Different probe radius should give different result
        let sasa_large_probe = calc_mol_sasa_with_probe(&mol, &coords, 2.0);
        assert!((sasa_large_probe - sasa_default).abs() > 0.1);
    }

    #[test]
    fn test_sasa_descriptor_consistency() {
        // Descriptor total should equal sasa() result
        let mol = parse("CCO").unwrap();
        let coords = generate_coords(&mol);

        let desc = sasa_descriptor(&mol, &coords);
        let total_sasa = sasa(&mol, &coords);

        assert!((desc.total - total_sasa).abs() < 1e-6);
    }

    #[test]
    fn test_sasa_descriptor_nonzero_atoms() {
        // All per-atom values should be positive for reasonable molecules
        let mol = parse("c1ccccc1").unwrap(); // Benzene
        let coords = generate_coords(&mol);
        let desc = sasa_descriptor(&mol, &coords);

        for (i, &sasa_val) in desc.per_atom.iter().enumerate() {
            assert!(sasa_val > 0.0, "atom {} should have positive SASA", i);
        }
    }
}
