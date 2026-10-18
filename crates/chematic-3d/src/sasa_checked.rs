//! Checked Shrake-Rupley SASA using caller-supplied coordinates and radii.

use crate::fmath::DetMath;
use chematic_core::Molecule;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SasaInputError {
    AtomLimit,
    Count {
        field: &'static str,
        expected: usize,
        actual: usize,
    },
    Coordinate {
        index: usize,
    },
    Radius {
        index: usize,
    },
    ProbeRadius,
    SpherePoints,
    Coincident {
        first: usize,
        second: usize,
    },
    WorkLimit,
    NumericalRange,
}

impl std::fmt::Display for SasaInputError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::AtomLimit => write!(f, "checked SASA supports at most 2048 atoms"),
            Self::Count {
                field,
                expected,
                actual,
            } => write!(
                f,
                "{field} must have one entry per atom: expected {expected}, got {actual}"
            ),
            Self::Coordinate { index } => write!(f, "non-finite coordinate at atom {index}"),
            Self::Radius { index } => {
                write!(f, "radius at atom {index} must be finite and positive")
            }
            Self::ProbeRadius => write!(f, "probe_radius must be finite and nonnegative"),
            Self::SpherePoints => write!(f, "sphere_points must be in 32..=65536"),
            Self::Coincident { first, second } => write!(
                f,
                "coincident equal inflated spheres at atoms {first} and {second}"
            ),
            Self::WorkLimit => write!(
                f,
                "checked SASA exceeds 50000000 potential point-neighbor checks"
            ),
            Self::NumericalRange => write!(
                f,
                "SASA area is outside representable finite positive float range"
            ),
        }
    }
}
impl std::error::Error for SasaInputError {}

#[derive(Debug, Clone, PartialEq)]
pub struct CheckedSasa {
    /// Compensated sum of per-atom accessible areas, in square Angstrom.
    pub total: f64,
    /// One area per input atom, in the original order.
    pub per_atom: Vec<f64>,
}

struct Sphere {
    center: [f64; 3],
    radius: f64,
    area: f64,
}
struct Neighbor {
    delta: [f64; 3],
    radius: f64,
    occluder: f64,
}
struct Occlusion {
    neighbors: Vec<Vec<Neighbor>>,
    covered: Vec<bool>,
}

fn validate(
    mol: &Molecule,
    coordinates: &[[f64; 3]],
    radii: &[f64],
    probe: f64,
    points: usize,
) -> Result<Vec<Sphere>, SasaInputError> {
    let n = mol.atom_count();
    if n > 2048 {
        return Err(SasaInputError::AtomLimit);
    }
    for (field, actual) in [("coordinates", coordinates.len()), ("radii", radii.len())] {
        if actual != n {
            return Err(SasaInputError::Count {
                field,
                expected: n,
                actual,
            });
        }
    }
    if !probe.is_finite() || probe < 0.0 {
        return Err(SasaInputError::ProbeRadius);
    }
    if !(32..=65536).contains(&points) {
        return Err(SasaInputError::SpherePoints);
    }
    let mut spheres = Vec::with_capacity(n);
    for (index, (&center, &radius)) in coordinates.iter().zip(radii).enumerate() {
        if center.iter().any(|v| !v.is_finite()) {
            return Err(SasaInputError::Coordinate { index });
        }
        if !radius.is_finite() || radius <= 0.0 {
            return Err(SasaInputError::Radius { index });
        }
        let radius = radius + probe;
        let area = 4.0 * std::f64::consts::PI * radius * radius;
        if !area.is_finite() || area <= 0.0 {
            return Err(SasaInputError::NumericalRange);
        }
        spheres.push(Sphere {
            center,
            radius,
            area,
        });
    }
    Ok(spheres)
}

fn neighbor(a: &Sphere, b: &Sphere) -> Neighbor {
    let delta = std::array::from_fn::<_, 3, _>(|k| a.center[k] - b.center[k]);
    let scale = delta
        .iter()
        .fold(a.radius.max(b.radius), |s, v| s.max(v.abs()));
    Neighbor {
        delta: delta.map(|v| v / scale),
        radius: a.radius / scale,
        occluder: b.radius / scale,
    }
}

fn prepare(spheres: &[Sphere], points: usize) -> Result<Occlusion, SasaInputError> {
    let n = spheres.len();
    let mut out = Occlusion {
        neighbors: (0..n).map(|_| Vec::new()).collect(),
        covered: vec![false; n],
    };
    let mut directed_pairs = 0usize;
    for i in 0..n {
        for j in i + 1..n {
            let (a, b) = (&spheres[i], &spheres[j]);
            let d = std::array::from_fn::<_, 3, _>(|k| a.center[k] - b.center[k]);
            let distance = libm::hypot(libm::hypot(d[0], d[1]), d[2]);
            if distance == 0.0 && a.radius == b.radius {
                return Err(SasaInputError::Coincident {
                    first: i,
                    second: j,
                });
            }
            if distance + a.radius <= b.radius {
                out.covered[i] = true;
                continue;
            }
            if distance + b.radius <= a.radius {
                out.covered[j] = true;
                continue;
            }
            if distance >= a.radius + b.radius {
                continue;
            }
            directed_pairs += 2;
            if directed_pairs
                .checked_mul(points)
                .is_none_or(|w| w > 50_000_000)
            {
                return Err(SasaInputError::WorkLimit);
            }
            out.neighbors[i].push(neighbor(a, b));
            out.neighbors[j].push(neighbor(b, a));
        }
    }
    Ok(out)
}

fn directions(points: usize) -> Vec<[f64; 3]> {
    let angle = std::f64::consts::PI * (3.0 - libm::sqrt(5.0));
    (0..points)
        .map(|i| {
            let y = 1.0 - 2.0 * (i as f64 + 0.5) / points as f64;
            let r = libm::sqrt((1.0 - y * y).max(0.0));
            let (s, c) = (angle * i as f64).dsin_cos();
            [c * r, y, s * r]
        })
        .collect()
}

fn exposed(point: &[f64; 3], neighbors: &[Neighbor]) -> bool {
    neighbors.iter().all(|n| {
        let d = std::array::from_fn::<_, 3, _>(|k| n.delta[k] + point[k] * n.radius);
        libm::hypot(libm::hypot(d[0], d[1]), d[2]) >= n.occluder
    })
}

/// Compute per-atom/total SASA with explicit Angstrom coordinates and positive
/// van der Waals radii, including explicit H atoms. Adds the nonnegative probe
/// radius exactly once. Does not generate coordinates or assign default radii.
///
/// Uses midpoint Fibonacci Shrake-Rupley sampling, not the FreeSASA backend.
/// Requires at most 2048 atoms, 32..=65536 points per sphere, and at most
/// 50 million potential point-neighbor checks. Identical coincident inflated
/// spheres are refused because per-atom surface ownership is ambiguous.
/// Nested unequal spheres are supported. All errors refuse partial output.
/// Empty molecules with empty input arrays return zero area.
///
/// Local normalized differences avoid absolute surface-point translation and
/// norm overflow/underflow; unrepresentable areas are errors. Sampling error
/// and orientation dependence remain, so this is not exact surface integration
/// or bitwise RDKit compatibility. Existing SASA entry points are unchanged.
pub fn sasa_with_radii_checked(
    mol: &Molecule,
    coordinates: &[[f64; 3]],
    radii: &[f64],
    probe_radius: f64,
    sphere_points: usize,
) -> Result<CheckedSasa, SasaInputError> {
    let spheres = validate(mol, coordinates, radii, probe_radius, sphere_points)?;
    let occlusion = prepare(&spheres, sphere_points)?;
    let directions = directions(sphere_points);
    let mut per_atom = Vec::with_capacity(spheres.len());
    let (mut total, mut compensation) = (0.0, 0.0);
    for (i, sphere) in spheres.iter().enumerate() {
        let count = if occlusion.covered[i] {
            0
        } else if occlusion.neighbors[i].is_empty() {
            sphere_points
        } else {
            directions
                .iter()
                .filter(|p| exposed(p, &occlusion.neighbors[i]))
                .count()
        };
        let area = sphere.area * (count as f64 / sphere_points as f64);
        if count > 0 && area == 0.0 {
            return Err(SasaInputError::NumericalRange);
        }
        let corrected = area - compensation;
        let next = total + corrected;
        if !next.is_finite() {
            return Err(SasaInputError::NumericalRange);
        }
        compensation = (next - total) - corrected;
        total = next;
        per_atom.push(area);
    }
    Ok(CheckedSasa { total, per_atom })
}

#[cfg(test)]
mod boundary_contract_tests {
    use super::*;
    use chematic_core::{Atom, Element, MoleculeBuilder};

    fn atoms(count: usize) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        for _ in 0..count {
            builder.add_atom(Atom::new(Element::C));
        }
        builder.build()
    }

    #[test]
    fn diagnostics_cover_every_input_contract() {
        let one = atoms(1);
        let cases = [
            sasa_with_radii_checked(&one, &[], &[1.0], 0.0, 64).unwrap_err(),
            sasa_with_radii_checked(&one, &[[0.0; 3]], &[], 0.0, 64).unwrap_err(),
            sasa_with_radii_checked(&one, &[[f64::NAN, 0.0, 0.0]], &[1.0], 0.0, 64).unwrap_err(),
            sasa_with_radii_checked(&one, &[[0.0; 3]], &[0.0], 0.0, 64).unwrap_err(),
            sasa_with_radii_checked(&one, &[[0.0; 3]], &[1.0], -0.1, 64).unwrap_err(),
            sasa_with_radii_checked(&one, &[[0.0; 3]], &[1.0], 0.0, 31).unwrap_err(),
            sasa_with_radii_checked(&one, &[[0.0; 3]], &[1.0e-200], 0.0, 64).unwrap_err(),
        ];
        for error in cases {
            assert!(!error.to_string().is_empty());
        }
        let too_many = atoms(2049);
        assert_eq!(
            sasa_with_radii_checked(&too_many, &vec![[0.0; 3]; 2049], &vec![1.0; 2049], 0.0, 64,)
                .unwrap_err(),
            SasaInputError::AtomLimit
        );
    }

    #[test]
    fn occlusion_handles_coincident_nested_disjoint_and_partial_spheres() {
        let two = atoms(2);
        assert!(matches!(
            sasa_with_radii_checked(&two, &[[0.0; 3]; 2], &[1.0; 2], 0.0, 64),
            Err(SasaInputError::Coincident {
                first: 0,
                second: 1
            })
        ));

        let nested = sasa_with_radii_checked(
            &two,
            &[[0.0, 0.0, 0.0], [0.25, 0.0, 0.0]],
            &[2.0, 0.5],
            0.0,
            128,
        )
        .unwrap();
        assert_eq!(nested.per_atom[1], 0.0);
        assert!(nested.per_atom[0] > 0.0);

        let disjoint = sasa_with_radii_checked(
            &two,
            &[[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]],
            &[1.0, 1.0],
            0.0,
            128,
        )
        .unwrap();
        assert!(disjoint.per_atom.iter().all(|area| *area > 12.0));

        let overlap = sasa_with_radii_checked(
            &two,
            &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            &[1.0, 1.0],
            0.0,
            128,
        )
        .unwrap();
        assert!(overlap.total < disjoint.total);
    }

    #[test]
    fn bounded_work_and_finite_sum_are_enforced() {
        let crowded = atoms(30);
        let coordinates: Vec<_> = (0..30).map(|i| [i as f64 * 0.01, 0.0, 0.0]).collect();
        assert_eq!(
            sasa_with_radii_checked(&crowded, &coordinates, &[10.0; 30], 0.0, 65_536).unwrap_err(),
            SasaInputError::WorkLimit
        );

        let large = atoms(20);
        let coordinates: Vec<_> = (0..20).map(|i| [i as f64 * 4.0e153, 0.0, 0.0]).collect();
        assert_eq!(
            sasa_with_radii_checked(&large, &coordinates, &[1.0e153; 20], 0.0, 32).unwrap_err(),
            SasaInputError::NumericalRange
        );
    }
}
