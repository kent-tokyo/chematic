// ---------------------------------------------------------------------------
// 7. TPSA
// ---------------------------------------------------------------------------

/// Compute the topological polar surface area (Å²) using the Ertl (2000) table.
///
/// Reference: P. Ertl, B. Rohde, P. Selzer, J. Med. Chem. 2000, 43, 3714-3717.
#[derive(Clone, Copy, Debug, Default)]
struct TpsaEnvironment {
    neighbors: u8,
    single: u8,
    double: u8,
    triple: u8,
    aromatic: u8,
    hydrogens: u8,
    charge: i8,
    in_three_membered_ring: bool,
}

fn tpsa_environment(mol: &Molecule, idx: AtomIdx, hydrogens: u8) -> TpsaEnvironment {
    let mut environment = TpsaEnvironment {
        hydrogens,
        charge: mol.atom(idx).charge,
        ..TpsaEnvironment::default()
    };
    let mut heavy_neighbors = smallvec::SmallVec::<[AtomIdx; 4]>::new();
    for (neighbor, bond_idx) in mol.neighbors(idx) {
        if mol.atom(neighbor).element.atomic_number() == 1 {
            continue;
        }
        environment.neighbors = environment.neighbors.saturating_add(1);
        heavy_neighbors.push(neighbor);
        match mol.bond(bond_idx).order {
            BondOrder::Single | BondOrder::Up | BondOrder::Down => {
                environment.single = environment.single.saturating_add(1);
            }
            BondOrder::Double => environment.double = environment.double.saturating_add(1),
            BondOrder::Triple => environment.triple = environment.triple.saturating_add(1),
            BondOrder::Aromatic => {
                environment.aromatic = environment.aromatic.saturating_add(1);
            }
            _ => {}
        }
    }
    environment.in_three_membered_ring = heavy_neighbors.iter().enumerate().any(|(i, &left)| {
        heavy_neighbors[i + 1..]
            .iter()
            .any(|&right| mol.bond_between(left, right).is_some())
    });
    environment
}

fn tpsa_nitrogen(environment: TpsaEnvironment) -> f64 {
    let e = environment;
    let bond_key = (
        e.neighbors,
        e.hydrogens,
        e.charge,
        e.single,
        e.double,
        e.triple,
        e.aromatic,
    );
    if bond_key == (2, 1, 0, 2, 0, 0, 0) {
        return if e.in_three_membered_ring {
            21.94
        } else {
            12.03
        };
    }
    if bond_key == (3, 0, 0, 3, 0, 0, 0) {
        return if e.in_three_membered_ring { 3.01 } else { 3.24 };
    }
    let key = (
        e.neighbors,
        e.hydrogens,
        e.charge,
        e.single,
        e.double,
        e.triple,
        e.aromatic,
        e.in_three_membered_ring,
    );
    let contribution = match key {
        (1, 0, 0, 0, 0, 1, 0, _) => Some(23.79),
        (1, 1, 0, 0, 1, 0, 0, _) => Some(23.85),
        (1, 2, 0, 1, 0, 0, 0, _) => Some(26.02),
        (1, 2, 1, 0, 1, 0, 0, _) => Some(25.59),
        (1, 3, 1, 1, 0, 0, 0, _) => Some(27.64),
        (2, 0, 0, 1, 1, 0, 0, _) => Some(12.36),
        (2, 0, 0, 0, 1, 1, 0, _) => Some(13.60),
        (2, 0, 1, 1, 0, 1, 0, _) => Some(4.36),
        (2, 1, 1, 1, 1, 0, 0, _) => Some(13.97),
        (2, 2, 1, 2, 0, 0, 0, _) => Some(16.61),
        (2, 0, 0, 0, 0, 0, 2, _) => Some(12.89),
        (2, 1, 0, 0, 0, 0, 2, _) => Some(15.79),
        (2, 1, 1, 0, 0, 0, 2, _) => Some(14.14),
        (3, 0, 0, 1, 2, 0, 0, _) => Some(11.68),
        (3, 0, 1, 2, 1, 0, 0, _) => Some(3.01),
        (3, 1, 1, 3, 0, 0, 0, _) => Some(4.44),
        (3, 0, 0, 0, 0, 0, 3, _) => Some(4.41),
        (3, 0, 0, 1, 0, 0, 2, _) => Some(4.93),
        (3, 0, 0, 0, 1, 0, 2, _) => Some(8.39),
        (3, 0, 1, 0, 0, 0, 3, _) => Some(4.10),
        (3, 0, 1, 1, 0, 0, 2, _) => Some(3.88),
        (4, 0, 1, 4, 0, 0, 0, _) => Some(0.0),
        _ => None,
    };
    contribution.unwrap_or_else(|| {
        (30.5 - f64::from(e.neighbors) * 8.2 + f64::from(e.hydrogens) * 1.5).max(0.0)
    })
}

fn tpsa_oxygen(environment: TpsaEnvironment) -> f64 {
    let e = environment;
    let contribution = match e.neighbors {
        1 if e.hydrogens == 0 && e.charge == 0 && e.double == 1 => Some(17.07),
        1 if e.hydrogens == 1 && e.charge == 0 && e.single == 1 => Some(20.23),
        1 if e.hydrogens == 0 && e.charge == -1 && e.single == 1 => Some(23.06),
        2 if e.hydrogens == 0 && e.charge == 0 && e.single == 2 && e.in_three_membered_ring => {
            Some(12.53)
        }
        2 if e.hydrogens == 0 && e.charge == 0 && e.single == 2 && !e.in_three_membered_ring => {
            Some(9.23)
        }
        2 if e.hydrogens == 0 && e.charge == 0 && e.aromatic == 2 => Some(13.14),
        _ => None,
    };
    contribution.unwrap_or_else(|| {
        (28.5 - f64::from(e.neighbors) * 8.6 + f64::from(e.hydrogens) * 1.5).max(0.0)
    })
}

fn tpsa_phosphorus(environment: TpsaEnvironment) -> f64 {
    let e = environment;
    match e.neighbors {
        2 if e.hydrogens == 0 && e.charge == 0 && e.single == 1 && e.double == 1 => 34.14,
        3 if e.hydrogens == 0 && e.charge == 0 && e.single == 3 => 13.59,
        3 if e.hydrogens == 1 && e.charge == 0 && e.single == 2 && e.double == 1 => 23.47,
        4 if e.hydrogens == 0 && e.charge == 0 && e.single == 3 && e.double == 1 => 9.81,
        _ => 0.0,
    }
}

fn tpsa_sulfur(environment: TpsaEnvironment) -> f64 {
    let e = environment;
    match e.neighbors {
        1 if e.hydrogens == 0 && e.charge == 0 && e.double == 1 => 32.09,
        1 if e.hydrogens == 1 && e.charge == 0 && e.single == 1 => 38.80,
        2 if e.hydrogens == 0 && e.charge == 0 && e.single == 2 => 25.30,
        2 if e.hydrogens == 0 && e.charge == 0 && e.aromatic == 2 => 28.24,
        3 if e.hydrogens == 0 && e.charge == 0 && e.aromatic == 2 && e.double == 1 => 21.70,
        3 if e.hydrogens == 0 && e.charge == 0 && e.single == 2 && e.double == 1 => 19.21,
        4 if e.hydrogens == 0 && e.charge == 0 && e.single == 2 && e.double == 2 => 8.38,
        _ => 0.0,
    }
}

/// Use the perceived (RDKit-parity) environment for a nitrogen that the input
/// spelled as non-aromatic but perception made aromatic. Imide nitrogens such
/// as phthalimide's are not aromatic under RDKit-parity perception, so they
/// keep the amide atom type without a special case.
fn use_perceived_nitrogen_environment(mol: &Molecule, mol_arom: &Molecule, idx: AtomIdx) -> bool {
    // RDKit types every nitrogen on the sanitized (aromaticity-perceived)
    // molecule. The earlier carbonyl/degree-2 restriction left Kekule inputs
    // such as N-alkyl carbazoles and pyridinium typed as aliphatic N; the
    // perceived environment matches RDKit on both aromatic and Kekule
    // spellings (ChEMBL 5k unchanged at 5000/5000, NCI 5k 4697 -> 4990).
    !mol.atom(idx).aromatic && mol_arom.atom(idx).aromatic
}

/// Aromaticity representation used by RDKit-compatible descriptor typing.
///
/// The general-purpose perception engine intentionally accepts a wider set of
/// conjugated ring systems. Descriptor atom types need RDKit's narrower
/// aromaticity boundary, especially for neutral divalent sulfur in fused
/// heterocycles. Keep the operation fail-safe for unusual inputs: a parity
/// perception failure falls back to the stable general model rather than
/// making descriptor evaluation fallible.
pub(crate) fn descriptor_aromaticity(mol: &Molecule) -> DescriptorView<'_> {
    // The RDKit-parity view is often an exact copy of `mol`; use `mol` (and
    // its caches) directly then.
    // The view also carries RDKit's sanitization rewrites (perchlorate,
    // organometallic dative bonds), so it is never `mol` when those apply.
    if chematic_perception::rdkit_parity_view_is_identity(mol) {
        return DescriptorView::Same(mol);
    }
    let shared = chematic_perception::apply_aromaticity_rdkit_parity_shared(mol);
    if let Ok(view) = shared.as_ref()
        && view
            .derived_if_computed::<Vec<bool>>(chematic_core::DerivedSlot::RdkitDativeBonds)
            .is_some_and(|flags| !flags.iter().any(|&dative| dative))
        && !view.bonds().any(|(_, bond)| {
            bond.order == chematic_core::BondOrder::Triple
                && view.atom(bond.atom1).aromatic
                && view.atom(bond.atom2).aromatic
        })
    {
        // This read-only view needs no descriptor-specific normalization,
        // including the no-further-dative marker normally seeded on a copy.
        // Externally seeded views without that marker keep the copy path.
        return DescriptorView::Shared(shared);
    }
    // Memoized on `mol`: TPSA, Crippen LogP/MR, HBA and friends all start from
    // the same perceived copy, and that copy carries its own SSSR cache.
    DescriptorView::Perceived(
        mol.derived(chematic_core::DerivedSlot::DescriptorAromatic, || {
            let mut view = chematic_perception::apply_aromaticity_rdkit_parity_experimental(mol)
                .unwrap_or_else(|_| {
                    let cleaned = chematic_perception::rdkit_sanitize_cleanup(mol);
                    chematic_perception::apply_aromaticity(cleaned.as_ref().unwrap_or(mol))
                });
            // A ring triple bond RDKit flags aromatic (`c1ccccc#1`) is kept
            // as a triple bond in the view; RDKit's descriptors read its
            // aromatic flag.
            let triples: Vec<chematic_core::BondIdx> = view
                .bonds()
                .filter(|(_, b)| {
                    b.order == chematic_core::BondOrder::Triple
                        && view.atom(b.atom1).aromatic
                        && view.atom(b.atom2).aromatic
                })
                .map(|(i, _)| i)
                .collect();
            for b in triples {
                view.set_bond_order(b, chematic_core::BondOrder::Aromatic);
            }
            view
        }),
    )
}

/// The descriptor aromatic view: `mol` itself or a perceived copy.
pub(crate) enum DescriptorView<'a> {
    Same(&'a Molecule),
    Perceived(std::sync::Arc<Molecule>),
    Shared(std::sync::Arc<Result<Molecule, chematic_perception::AromaticityError>>),
}

impl std::ops::Deref for DescriptorView<'_> {
    type Target = Molecule;
    fn deref(&self) -> &Molecule {
        match self {
            DescriptorView::Same(m) => m,
            DescriptorView::Perceived(m) => m,
            // safe: `Shared` is only constructed after checking `Ok` above;
            // the immutable Arc result cannot change after that check.
            DescriptorView::Shared(m) => m.as_ref().as_ref().expect("validated aromatic view"),
        }
    }
}

fn tpsa_contributions(mol: &Molecule) -> Vec<f64> {
    let mol_arom = descriptor_aromaticity(mol);
    let mut contributions = vec![0.0; mol.atom_count()];
    for (idx, original_atom) in mol.atoms() {
        let atomic_number = original_atom.element.atomic_number();
        contributions[idx.0 as usize] = match atomic_number {
            7 => {
                let hydrogens = descriptor_attached_hcount(mol, idx);
                let environment_mol = if use_perceived_nitrogen_environment(mol, &mol_arom, idx) {
                    &mol_arom
                } else {
                    mol
                };
                tpsa_nitrogen(tpsa_environment(environment_mol, idx, hydrogens))
            }
            8 => {
                let hydrogens = descriptor_attached_hcount(&mol_arom, idx);
                tpsa_oxygen(tpsa_environment(&mol_arom, idx, hydrogens))
            }
            15 => {
                let hydrogens = descriptor_attached_hcount(&mol_arom, idx);
                tpsa_phosphorus(tpsa_environment(&mol_arom, idx, hydrogens))
            }
            16 => {
                let hydrogens = descriptor_attached_hcount(&mol_arom, idx);
                tpsa_sulfur(tpsa_environment(&mol_arom, idx, hydrogens))
            }
            _ => 0.0,
        };
    }
    contributions
}

/// Topological Polar Surface Area (Ertl 2000).
///
/// Sum of Ertl atom-type contributions for N, O, S, and P atoms.
/// Hydrogen atoms are implicit (not computed). Values match RDKit defaults with calibrations
/// for secondary amide N (12.03 Å²), aromatic N (15.79 Å²), and phosphorus atoms.
///
/// **Note on S/P:** S and P contributions are **included by default** (equivalent to
/// RDKit `Descriptors.TPSA(mol, includeSandP=True)`). RDKit's default
/// `Descriptors.TPSA(mol)` excludes S and P; results will differ for molecules
/// containing these elements.
pub fn tpsa(mol: &Molecule) -> f64 {
    mol.derived_scalar(chematic_core::derived_cache::DerivedScalar::Tpsa, || {
        tpsa_contributions(mol).into_iter().sum::<f64>()
    })
}

/// RDKit-default TPSA: `rdMolDescriptors.CalcTPSA(mol)` / `Descriptors.TPSA(mol)`
/// with `includeSandP=False`, i.e. the N and O contributions of [`tpsa`] only.
pub fn rdkit_tpsa(mol: &Molecule) -> f64 {
    chematic_smiles::with_rdkit_model_scalar(
        mol,
        chematic_core::derived_cache::DerivedScalar::RdkitTpsa,
        rdkit_tpsa_impl,
    )
}

fn rdkit_tpsa_impl(mol: &Molecule) -> f64 {
    // This default model sums only N/O contributions. With neither element,
    // the original fold adds only positive zeros, without needing perception.
    if !mol
        .atoms()
        .any(|(_, atom)| matches!(atom.element.atomic_number(), 7 | 8))
    {
        return 0.0;
    }
    let mol_arom = descriptor_aromaticity(mol);
    // RDKit's default excludes S/P. Compute only N/O instead of allocating a
    // full contribution vector and typing atoms whose values are discarded.
    // Iteration and accumulation order remain atom-index order, exactly as in
    // `tpsa_contributions`, so floating-point results are unchanged.
    mol.atoms().fold(0.0, |acc, (idx, atom)| {
        let contribution = match atom.element.atomic_number() {
            7 => {
                let hydrogens = descriptor_attached_hcount(mol, idx);
                let environment_mol = if use_perceived_nitrogen_environment(mol, &mol_arom, idx) {
                    &mol_arom
                } else {
                    mol
                };
                tpsa_nitrogen(tpsa_environment(environment_mol, idx, hydrogens))
            }
            8 => {
                let hydrogens = descriptor_attached_hcount(&mol_arom, idx);
                tpsa_oxygen(tpsa_environment(&mol_arom, idx, hydrogens))
            }
            _ => 0.0,
        };
        acc + contribution
    })
}

#[cfg(test)]
#[path = "descriptor_view_tests.rs"]
mod descriptor_view_tests;

#[cfg(test)]
#[path = "tpsa_zero_tests.rs"]
mod tpsa_zero_tests;
