//! Electronegativity-equalization-method (EEM) partial charges.
//!
//! This implements the constrained linear system used by RDKit's
//! `CalcEEMcharges`, with the charge-conservation row oriented correctly.
//! RDKit issue #9688 documents that its current matrix mapping reverses the
//! requested molecular charge for ions.  The implementation here uses an
//! explicit row-major solver, so `sum(q_i)` is the molecule's formal charge.

use chematic_core::{AtomIdx, BondOrder, Coords3D, Molecule, apply_kekule, kekulize};

const KAPPA: f64 = 0.5125;
/// The dense solve is O(n^3); keep hostile inputs from monopolising a worker.
pub const MAX_EEM_ATOMS: usize = 512;

// CCD_gen_DE_RMSD_B3LYP_6311G_NPA parameters published with the NEEMP
// reference model (https://doi.org/10.1186/s13321-016-0171-1) and used by
// RDKit. Atomic numbers not listed by that model are intentionally unsupported
// rather than silently receiving zero charge.
const A1: [f64; 36] = [
    0.0, 2.5473, 0.0, 0.0, 0.0, 0.0, 2.7221, 2.9750, 3.1503, 2.9976, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    2.6511, 2.7026, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 2.6263,
];
const A2: [f64; 36] = [
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.7667, 2.8895, 3.0486, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.2933,
    2.6471, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0,
];
const A3: [f64; 36] = [
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.6944, 3.0240, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
];
const B1: [f64; 36] = [
    0.0, 1.1641, 0.0, 0.0, 0.0, 0.0, 0.6403, 0.9083, 1.0577, 0.9983, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.4897, 1.1537, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.1105,
];
const B2: [f64; 36] = [
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6513, 0.6647, 0.8410, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5759,
    0.4512, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0,
];
const B3: [f64; 36] = [
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6776, 1.4240, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
];

#[derive(Debug, Clone, PartialEq)]
pub enum EemError {
    EmptyMolecule,
    TooManyAtoms {
        count: usize,
        limit: usize,
    },
    CoordinateCount {
        atoms: usize,
        coordinates: usize,
    },
    NonFiniteCoordinate {
        atom: usize,
    },
    CoincidentAtoms {
        first: usize,
        second: usize,
    },
    UnsupportedAtom {
        atom: usize,
        atomic_number: u8,
        bond_type: u8,
    },
    Kekulization(String),
    SingularSystem,
}

impl core::fmt::Display for EemError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::EmptyMolecule => write!(f, "EEM requires at least one atom"),
            Self::TooManyAtoms { count, limit } => {
                write!(f, "EEM atom count {count} exceeds limit {limit}")
            }
            Self::CoordinateCount { atoms, coordinates } => {
                write!(f, "EEM expected {atoms} coordinates, got {coordinates}")
            }
            Self::NonFiniteCoordinate { atom } => {
                write!(f, "EEM coordinate for atom {atom} is not finite")
            }
            Self::CoincidentAtoms { first, second } => {
                write!(
                    f,
                    "EEM atoms {first} and {second} have coincident coordinates"
                )
            }
            Self::UnsupportedAtom {
                atom,
                atomic_number,
                bond_type,
            } => write!(
                f,
                "EEM has no parameters for atom {atom} (atomic number {atomic_number}, bond type {bond_type})"
            ),
            Self::Kekulization(detail) => write!(f, "EEM kekulization failed: {detail}"),
            Self::SingularSystem => write!(f, "EEM linear system is singular"),
        }
    }
}

impl std::error::Error for EemError {}

fn atom_type(mol: &Molecule, atom: AtomIdx) -> u8 {
    mol.neighbors(atom)
        .map(|(_, bond)| match mol.bond(bond).order {
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            BondOrder::Quadruple => 4,
            _ => 1,
        })
        .max()
        .unwrap_or(1)
}

fn parameters(atomic_number: u8, bond_type: u8) -> Option<(f64, f64)> {
    let index = atomic_number as usize;
    if index >= A1.len() {
        return None;
    }
    let (a, b) = match bond_type {
        1 => (A1[index], B1[index]),
        2 => (A2[index], B2[index]),
        3 => (A3[index], B3[index]),
        _ => return None,
    };
    (a != 0.0 && b != 0.0).then_some((a, b))
}

/// Solve a dense linear system in place with partial pivoting.
fn solve(mut matrix: Vec<Vec<f64>>, mut rhs: Vec<f64>) -> Option<Vec<f64>> {
    let n = rhs.len();
    for column in 0..n {
        let pivot = (column..n)
            .max_by(|&a, &b| matrix[a][column].abs().total_cmp(&matrix[b][column].abs()))?;
        if matrix[pivot][column].abs() < 1e-12 {
            return None;
        }
        matrix.swap(column, pivot);
        rhs.swap(column, pivot);
        let diagonal = matrix[column][column];
        for row in (column + 1)..n {
            let factor = matrix[row][column] / diagonal;
            if factor == 0.0 {
                continue;
            }
            matrix[row][column] = 0.0;
            for k in (column + 1)..n {
                matrix[row][k] -= factor * matrix[column][k];
            }
            rhs[row] -= factor * rhs[column];
        }
    }
    let mut answer = vec![0.0; n];
    for row in (0..n).rev() {
        let known: f64 = ((row + 1)..n)
            .map(|column| matrix[row][column] * answer[column])
            .sum();
        answer[row] = (rhs[row] - known) / matrix[row][row];
    }
    answer
        .iter()
        .all(|value| value.is_finite())
        .then_some(answer)
}

/// Compute per-atom EEM partial charges from one 3D conformer.
///
/// The returned values conserve formal charge, including for ions. Supported
/// elements are H, C, N, O, F, P, S, Cl and Br within the published model's
/// available single/double/triple-bond atom types.
pub fn eem_charges(mol: &Molecule, coords: &Coords3D) -> Result<Vec<f64>, EemError> {
    let n = mol.atom_count();
    if n == 0 {
        return Err(EemError::EmptyMolecule);
    }
    if n > MAX_EEM_ATOMS {
        return Err(EemError::TooManyAtoms {
            count: n,
            limit: MAX_EEM_ATOMS,
        });
    }
    if coords.atom_count() != n {
        return Err(EemError::CoordinateCount {
            atoms: n,
            coordinates: coords.atom_count(),
        });
    }
    for (atom, point) in coords.points.iter().enumerate() {
        if !point.is_finite() {
            return Err(EemError::NonFiniteCoordinate { atom });
        }
    }

    let kekule = kekulize(mol).map_err(|error| EemError::Kekulization(error.detail))?;
    let view = apply_kekule(mol, &kekule);
    let size = n + 1;
    let mut matrix = vec![vec![0.0; size]; size];
    let mut rhs = vec![0.0; size];

    for i in 0..n {
        let atom = view.atom(AtomIdx(i as u32));
        let bond_type = atom_type(&view, AtomIdx(i as u32));
        let (electronegativity, hardness) = parameters(atom.element.atomic_number(), bond_type)
            .ok_or(EemError::UnsupportedAtom {
                atom: i,
                atomic_number: atom.element.atomic_number(),
                bond_type,
            })?;
        matrix[i][i] = hardness;
        rhs[i] = -electronegativity;
        matrix[i][n] = 1.0;
        matrix[n][i] = 1.0;
        for j in (i + 1)..n {
            let distance = coords.points[i].distance(&coords.points[j]);
            if distance <= 1e-12 {
                return Err(EemError::CoincidentAtoms {
                    first: i,
                    second: j,
                });
            }
            let interaction = KAPPA / distance;
            matrix[i][j] = interaction;
            matrix[j][i] = interaction;
        }
    }
    let formal_charge = view.atoms().map(|(_, atom)| f64::from(atom.charge)).sum();
    rhs[n] = formal_charge;

    let solution = solve(matrix, rhs).ok_or(EemError::SingularSystem)?;
    let mut charges = solution[..n].to_vec();
    // Remove the last few ulps of solver drift without changing relative charges.
    let residual = formal_charge - charges.iter().sum::<f64>();
    charges[n - 1] += residual;
    Ok(charges)
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_core::{Atom, Element, MoleculeBuilder, Point3};
    use chematic_smiles::parse;

    fn coords(rows: &[[f64; 3]]) -> Coords3D {
        Coords3D {
            points: rows
                .iter()
                .map(|row| Point3::new(row[0], row[1], row[2]))
                .collect(),
        }
    }

    #[test]
    fn charged_molecules_conserve_the_requested_charge() {
        let hydroxide = parse("[O-][H]").unwrap();
        let q = eem_charges(&hydroxide, &coords(&[[0.0, 0.0, 0.0], [0.97, 0.0, 0.0]])).unwrap();
        assert!((q.iter().sum::<f64>() + 1.0).abs() < 1e-12, "{q:?}");
        assert!(q[0] < -1.0 && q[1] > 0.0, "{q:?}");

        let ammonium = parse("[N+]([H])([H])([H])[H]").unwrap();
        let a = 1.01 / 3.0_f64.sqrt();
        let q = eem_charges(
            &ammonium,
            &coords(&[
                [0.0, 0.0, 0.0],
                [a, a, a],
                [a, -a, -a],
                [-a, a, -a],
                [-a, -a, a],
            ]),
        )
        .unwrap();
        assert!((q.iter().sum::<f64>() - 1.0).abs() < 1e-12, "{q:?}");
        assert!(
            q[0] < 0.0 && q[1..].iter().all(|value| *value > 0.0),
            "{q:?}"
        );
        assert!(
            q[1..]
                .windows(2)
                .all(|pair| (pair[0] - pair[1]).abs() < 1e-12)
        );
    }

    #[test]
    fn neutral_reference_matches_rdkit_parameterization() {
        let water = parse("O([H])[H]").unwrap();
        let q = eem_charges(
            &water,
            &coords(&[[0.0, 0.0, 0.0], [0.9572, 0.0, 0.0], [-0.239, 0.927, 0.0]]),
        )
        .unwrap();
        assert!(q.iter().sum::<f64>().abs() < 1e-12, "{q:?}");
        let expected = [-0.8167316014186569, 0.40839732995340194, 0.4083342714652547];
        assert!(
            q.iter()
                .zip(expected)
                .all(|(actual, expected)| (actual - expected).abs() < 1e-12),
            "{q:?}"
        );
    }

    #[test]
    fn rejects_unsupported_elements_and_degenerate_geometry() {
        let iodine = parse("I").unwrap();
        assert!(matches!(
            eem_charges(&iodine, &coords(&[[0.0, 0.0, 0.0]])),
            Err(EemError::UnsupportedAtom {
                atomic_number: 53,
                ..
            })
        ));
        let quadruple = parse("[C]$[C]").unwrap();
        assert!(matches!(
            eem_charges(&quadruple, &coords(&[[0.0, 0.0, 0.0], [1.2, 0.0, 0.0]])),
            Err(EemError::UnsupportedAtom { bond_type: 4, .. })
        ));
        let hydrogen = parse("[H][H]").unwrap();
        assert!(matches!(
            eem_charges(&hydrogen, &coords(&[[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])),
            Err(EemError::CoincidentAtoms { .. })
        ));
    }

    #[test]
    fn rejects_bounded_and_malformed_inputs_with_specific_errors() {
        let empty = MoleculeBuilder::new().build();
        assert_eq!(
            eem_charges(&empty, &coords(&[])),
            Err(EemError::EmptyMolecule)
        );

        let mut builder = MoleculeBuilder::new();
        for _ in 0..=MAX_EEM_ATOMS {
            builder.add_atom(Atom::new(Element::H));
        }
        assert_eq!(
            eem_charges(&builder.build(), &coords(&[])),
            Err(EemError::TooManyAtoms {
                count: MAX_EEM_ATOMS + 1,
                limit: MAX_EEM_ATOMS,
            })
        );

        let hydrogen = parse("[H][H]").unwrap();
        assert_eq!(
            eem_charges(&hydrogen, &coords(&[[0.0, 0.0, 0.0]])),
            Err(EemError::CoordinateCount {
                atoms: 2,
                coordinates: 1,
            })
        );
        assert_eq!(
            eem_charges(&hydrogen, &coords(&[[0.0, 0.0, 0.0], [f64::NAN, 0.0, 0.0]])),
            Err(EemError::NonFiniteCoordinate { atom: 1 })
        );

        let singular_distance = KAPPA / B1[1];
        assert_eq!(
            eem_charges(
                &hydrogen,
                &coords(&[[0.0, 0.0, 0.0], [singular_distance, 0.0, 0.0]])
            ),
            Err(EemError::SingularSystem)
        );
    }

    #[test]
    fn covers_supported_multiple_bond_types_and_kekulization_failure() {
        for smiles in ["C=C", "C#N"] {
            let mol = parse(smiles).unwrap();
            let charges = eem_charges(&mol, &coords(&[[0.0, 0.0, 0.0], [1.3, 0.0, 0.0]])).unwrap();
            assert!(charges.iter().sum::<f64>().abs() < 1e-12, "{smiles}");
        }

        let mut builder = MoleculeBuilder::new();
        let atoms = [
            builder.add_atom(Atom::aromatic(Element::C)),
            builder.add_atom(Atom::aromatic(Element::C)),
            builder.add_atom(Atom::aromatic(Element::C)),
        ];
        for i in 0..3 {
            builder
                .add_bond(atoms[i], atoms[(i + 1) % 3], BondOrder::Aromatic)
                .unwrap();
        }
        assert!(matches!(
            eem_charges(
                &builder.build(),
                &coords(&[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 0.8, 0.0]])
            ),
            Err(EemError::Kekulization(_))
        ));
    }

    #[test]
    fn error_messages_and_solver_failure_paths_are_stable() {
        let errors = [
            EemError::EmptyMolecule,
            EemError::TooManyAtoms { count: 2, limit: 1 },
            EemError::CoordinateCount {
                atoms: 2,
                coordinates: 1,
            },
            EemError::NonFiniteCoordinate { atom: 1 },
            EemError::CoincidentAtoms {
                first: 0,
                second: 1,
            },
            EemError::UnsupportedAtom {
                atom: 0,
                atomic_number: 53,
                bond_type: 1,
            },
            EemError::Kekulization("fixture".into()),
            EemError::SingularSystem,
        ];
        for error in errors {
            assert!(error.to_string().starts_with("EEM"));
        }

        assert_eq!(
            solve(vec![vec![1.0, 0.0], vec![0.0, 1.0]], vec![2.0, 3.0]),
            Some(vec![2.0, 3.0])
        );
        assert_eq!(solve(vec![vec![0.0]], vec![1.0]), None);
        assert_eq!(solve(vec![vec![f64::NAN]], vec![1.0]), None);
    }
}
