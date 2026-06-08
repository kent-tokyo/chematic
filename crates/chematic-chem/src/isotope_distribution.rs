//! Isotope distribution (isotope envelope) calculation.
//!
//! [`isotope_distribution`] computes the theoretical mass spectrum of a
//! molecule by convolving each atom's natural isotope abundances.  The result
//! is a list of `(m/z, relative_intensity)` pairs sorted by mass, normalised
//! so the most abundant peak has intensity 1.0.
//!
//! The `resolution` parameter merges peaks whose masses differ by less than
//! `resolution` Da — use `0.0` to keep every peak separate, `0.01` for
//! nominal-mass spectra, or `~0.001` for high-resolution instruments.

use chematic_core::{Element, Molecule, implicit_hcount};

// ---------------------------------------------------------------------------
// Isotope tables
// ---------------------------------------------------------------------------

/// `(mass_Da, natural_abundance_fraction)` for each isotope of an element.
/// Only isotopes with abundance > 0.0001 (0.01 %) are listed.
fn isotopes(element: Element) -> &'static [(f64, f64)] {
    match element.atomic_number() {
        1 => &[(1.00783, 0.99985), (2.01410, 0.00015)], // H
        6 => &[(12.0000, 0.9893), (13.00335, 0.0107)],  // C
        7 => &[(14.0031, 0.99636), (15.0001, 0.00364)], // N
        8 => &[(15.9949, 0.99757), (16.9991, 0.00038), (17.9992, 0.00205)], // O
        9 => &[(18.9984, 1.0)],                         // F
        11 => &[(22.9898, 1.0)],                        // Na
        14 => &[(27.9769, 0.9223), (28.9765, 0.0467), (29.9738, 0.0310)], // Si
        15 => &[(30.9738, 1.0)],                        // P
        16 => &[
            (31.9721, 0.9493),
            (32.9715, 0.0076),
            (33.9679, 0.0429),
            (35.9671, 0.0002),
        ], // S
        17 => &[(34.9689, 0.7576), (36.9659, 0.2424)],  // Cl
        19 => &[(38.9637, 0.9326), (40.9618, 0.0673)],  // K (approx)
        33 => &[(74.9216, 1.0)],                        // As
        34 => &[
            (73.9225, 0.0089),
            (75.9192, 0.0937),
            (76.9199, 0.0763),
            (77.9173, 0.2377),
            (79.9165, 0.4961),
            (81.9167, 0.0873),
        ], // Se
        35 => &[(78.9183, 0.5069), (80.9163, 0.4931)],  // Br
        53 => &[(126.9045, 1.0)],                       // I
        // Default: treat as monoisotopic (100 % at the nominal atomic number mass).
        n => {
            // Static fallback — only the most common elements reach here in practice.
            // We return a single-isotope slice stored as a leaked box to satisfy
            // the lifetime requirement.  This path is rarely hit.
            let _ = n; // suppress unused warning
            &[(0.0, 1.0)] // placeholder; overridden below
        }
    }
}

/// Return the monoisotopic mass of an element for the fallback case.
fn mono_mass_fallback(an: u8) -> f64 {
    match an {
        1 => 1.00783,
        6 => 12.0000,
        7 => 14.0031,
        8 => 15.9949,
        9 => 18.9984,
        14 => 27.9769,
        15 => 30.9738,
        16 => 31.9721,
        17 => 34.9689,
        35 => 78.9183,
        34 => 79.9165,
        53 => 126.9045,
        n => n as f64,
    }
}

/// Retrieve isotope distribution for an element, handling the fallback case.
fn isotope_list(element: Element) -> Vec<(f64, f64)> {
    let an = element.atomic_number();
    let table = isotopes(element);
    // The fallback branch returns `(0.0, 1.0)` — detect and replace.
    if table.len() == 1 && table[0].0 == 0.0 {
        vec![(mono_mass_fallback(an), 1.0)]
    } else {
        table.to_vec()
    }
}

// ---------------------------------------------------------------------------
// Convolution
// ---------------------------------------------------------------------------

/// Convolve two peak lists: `(mass, intensity)` pairs.
///
/// The resulting list has one entry for every combination of (m_a, m_b):
/// mass = m_a + m_b, intensity = i_a × i_b.
fn convolve(a: &[(f64, f64)], b: &[(f64, f64)]) -> Vec<(f64, f64)> {
    let mut out = Vec::with_capacity(a.len() * b.len());
    for &(ma, ia) in a {
        for &(mb, ib) in b {
            out.push((ma + mb, ia * ib));
        }
    }
    out
}

/// Merge peaks whose masses differ by less than `resolution` Da, summing
/// their intensities.  The merged peak mass is the intensity-weighted centroid.
fn merge_peaks(mut peaks: Vec<(f64, f64)>, resolution: f64) -> Vec<(f64, f64)> {
    if peaks.is_empty() {
        return peaks;
    }
    peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    if resolution <= 0.0 {
        return peaks;
    }

    let mut merged: Vec<(f64, f64)> = Vec::new();
    let mut cur_mass = peaks[0].0;
    let mut cur_int = peaks[0].1;

    for &(m, i) in &peaks[1..] {
        if m - cur_mass < resolution {
            // Intensity-weighted centroid merge.
            cur_mass = (cur_mass * cur_int + m * i) / (cur_int + i);
            cur_int += i;
        } else {
            merged.push((cur_mass, cur_int));
            cur_mass = m;
            cur_int = i;
        }
    }
    merged.push((cur_mass, cur_int));
    merged
}

/// Normalise peak intensities so the most abundant peak has intensity 1.0.
fn normalise(mut peaks: Vec<(f64, f64)>) -> Vec<(f64, f64)> {
    let max_int = peaks.iter().map(|&(_, i)| i).fold(0.0_f64, f64::max);
    if max_int > 0.0 {
        for (_, i) in &mut peaks {
            *i /= max_int;
        }
    }
    peaks
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Compute the theoretical isotope distribution of `mol`.
///
/// Returns a list of `(m/z, relative_intensity)` pairs sorted by mass.
/// The most abundant peak has intensity 1.0.
///
/// `resolution` (Da) controls peak merging:
/// - `0.0`  — keep every isotopologue as a separate entry.
/// - `0.01` — merge peaks within 0.01 Da (typical for low-res spectra).
/// - `1.0`  — nominal-mass / unit-resolution spectrum.
///
/// Explicit isotope labels on atoms (`atom.isotope: Some(mass_number)`) are
/// treated as pure single-isotope species at integer mass (same convention as
/// [`exact_mass`]).
pub fn isotope_distribution(mol: &Molecule, resolution: f64) -> Vec<(f64, f64)> {
    // Start with a single peak at mass 0.
    let mut dist: Vec<(f64, f64)> = vec![(0.0, 1.0)];

    for (idx, atom) in mol.atoms() {
        // Atom's own isotope distribution.
        let atom_dist: Vec<(f64, f64)> = match atom.isotope {
            Some(iso) => vec![(iso as f64, 1.0)], // explicit label → monoisotopic
            None => isotope_list(atom.element),
        };
        dist = convolve(&dist, &atom_dist);

        // Implicit hydrogens.
        let h_count = implicit_hcount(mol, idx) as usize;
        if h_count > 0 {
            let h_dist = isotope_list(chematic_core::Element::H);
            for _ in 0..h_count {
                dist = convolve(&dist, &h_dist);
            }
        }
    }

    // Prune negligible peaks (< 0.0001 % of total intensity) to keep the list manageable.
    let total: f64 = dist.iter().map(|&(_, i)| i).sum();
    if total > 0.0 {
        dist.retain(|(_, i)| *i / total >= 1e-6);
    }

    let dist = merge_peaks(dist, resolution);
    normalise(dist)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    /// Find the peak with the highest intensity.
    fn base_peak(dist: &[(f64, f64)]) -> (f64, f64) {
        dist.iter()
            .cloned()
            .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
            .unwrap()
    }

    #[test]
    fn test_water_monoisotopic_dominant() {
        // H2O: monoisotopic peak at 18.010 Da should be the largest.
        let mol = parse("O").unwrap();
        let dist = isotope_distribution(&mol, 0.0);
        assert!(!dist.is_empty());
        let (mass, _) = base_peak(&dist);
        assert!(
            (mass - 18.0106).abs() < 0.01,
            "H2O base peak ≈ 18.011 Da, got {mass:.4}"
        );
    }

    #[test]
    fn test_chloromethane_two_major_peaks() {
        // CH3Cl: ³⁵Cl (75.76 %) and ³⁷Cl (24.24 %) give two peaks separated by 2 Da.
        let mol = parse("CCl").unwrap();
        let dist = isotope_distribution(&mol, 0.01);
        // At least two peaks with significant intensity.
        let significant: Vec<_> = dist.iter().filter(|&&(_, i)| i > 0.1).collect();
        assert!(
            significant.len() >= 2,
            "CH3Cl should show M and M+2 peaks from Cl isotopes"
        );
    }

    #[test]
    fn test_bromobenzene_roughly_equal_peaks() {
        // C6H5Br: ⁷⁹Br / ⁸¹Br ≈ 50:50.  M and M+2 peaks should be close.
        let mol = parse("c1ccccc1Br").unwrap();
        let dist = isotope_distribution(&mol, 0.01);
        let top2: Vec<_> = {
            let mut sorted = dist.clone();
            sorted.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
            sorted.into_iter().take(2).collect()
        };
        assert_eq!(top2.len(), 2);
        let ratio = top2[0].1 / top2[1].1;
        assert!(
            ratio < 1.25,
            "Br isotopes should give near-equal peaks; ratio was {ratio:.3}"
        );
    }

    #[test]
    fn test_benzene_monoisotopic() {
        // C6H6: monoisotopic at 78.047 Da.
        let mol = parse("c1ccccc1").unwrap();
        let dist = isotope_distribution(&mol, 0.0);
        let (mass, _) = base_peak(&dist);
        assert!(
            (mass - 78.047).abs() < 0.01,
            "C6H6 M+ ≈ 78.047, got {mass:.4}"
        );
    }

    #[test]
    fn test_normalised_base_peak_is_one() {
        let mol = parse("CC(=O)O").unwrap(); // acetic acid
        let dist = isotope_distribution(&mol, 0.01);
        let max_int = dist.iter().map(|&(_, i)| i).fold(0.0_f64, f64::max);
        assert!(
            (max_int - 1.0).abs() < 1e-9,
            "base peak should be normalised to 1.0"
        );
    }

    #[test]
    fn test_resolution_merging() {
        let mol = parse("O").unwrap();
        let dist_fine = isotope_distribution(&mol, 0.0);
        let dist_coarse = isotope_distribution(&mol, 1.0);
        assert!(
            dist_coarse.len() <= dist_fine.len(),
            "coarser resolution should produce fewer or equal peaks"
        );
    }
}
