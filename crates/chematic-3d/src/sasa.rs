//! Solvent-Accessible Surface Area (SASA) calculation using Shrake-Rupley algorithm.

use chematic_core::{AtomIdx, Molecule};

use crate::coords::{Coords3D, Point3};

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
        let exposed_count = count_exposed_points(
            mol,
            coords,
            pos_i,
            idx,
            probe_radius,
            sphere_points,
        );

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

        let exposed_count = count_exposed_points(
            mol,
            coords,
            pos_i,
            idx,
            probe_radius,
            sphere_points,
        );

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
        return vec![center];  // Avoid division by zero when num_points < 2
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
        assert!((sum - total).abs() < 1e-6, "per-atom sum should match total");
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
}
