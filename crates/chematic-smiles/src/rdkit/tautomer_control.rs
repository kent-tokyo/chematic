//! Checked enumeration limits and synchronous cancellation checkpoints.
use super::{
    Molecule, RdkitSmilesError, RdkitSmilesParams, RdkitTautomerStatus, rdkit_mol_from_smiles,
    tautomer, write,
};

/// The result of [`rdkit_enumerate_tautomers`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RdkitTautomerEnumeration {
    /// `Chem.MolToSmiles(t)` of every tautomer, in the result's order (its
    /// map keys' order).
    pub smiles: Vec<String>,
    /// The result's map keys (`TautomerEnumeratorResult.smiles`).
    pub keys: Vec<String>,
    /// `TautomerEnumeratorResult.status`.
    pub status: RdkitTautomerStatus,
}

/// `rdMolStandardize.TautomerEnumerator().Enumerate(m)` for
/// `m = Chem.MolFromSmiles(s)` (RDKit 2026.03.1, default
/// `CleanupParameters`): the tautomers' SMILES and the result status.
pub fn rdkit_enumerate_tautomers(
    mol: &Molecule,
) -> Result<RdkitTautomerEnumeration, RdkitSmilesError> {
    let (m, _) = rdkit_mol_from_smiles(mol)?;
    let res = tautomer::enumerate(&m, &tautomer::Settings::default())?;
    let mut smiles = Vec::with_capacity(res.tautomers.len());
    for t in res.tautomers.values() {
        smiles.push(write::mol_to_smiles(&t.mol, &RdkitSmilesParams::default())?);
    }
    Ok(RdkitTautomerEnumeration {
        smiles,
        keys: res.tautomers.keys().cloned().collect(),
        status: res.status,
    })
}

/// Immutable snapshot before a matched transform is applied. Keys are the
/// current enumeration map keys, not freshly canonicalized final SMILES.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RdkitTautomerProgress {
    /// Number of transform rules with a nonempty match set so far.
    pub transforms: usize,
    /// Current unique tautomer keys, in lexicographic order.
    pub keys: Vec<String>,
}

/// Enumerate with positive limits and a synchronous callback before each
/// matched transform. Return false to stop with `Canceled` and a normalized
/// partial result. Limits take precedence over callbacks. No checkpoint is
/// emitted if no transform matches. Panics from a Rust callback propagate.
pub fn rdkit_enumerate_tautomers_with_callback(
    mol: &Molecule,
    max_tautomers: usize,
    max_transforms: usize,
    callback: &mut dyn FnMut(&RdkitTautomerProgress) -> bool,
) -> Result<RdkitTautomerEnumeration, RdkitSmilesError> {
    if max_tautomers == 0 || max_transforms == 0 {
        return Err(RdkitSmilesError::Unsupported(
            "tautomer limits must be positive".into(),
        ));
    }
    let (m, _) = rdkit_mol_from_smiles(mol)?;
    let settings = tautomer::Settings {
        max_tautomers,
        max_transforms,
        ..Default::default()
    };
    let mut checkpoint = |transforms, result: &tautomer::EnumerateResult| {
        callback(&RdkitTautomerProgress {
            transforms,
            keys: result.tautomers.keys().cloned().collect(),
        })
    };
    let res = tautomer::enumerate_with_callback(&m, &settings, Some(&mut checkpoint))?;
    let smiles = res
        .tautomers
        .values()
        .map(|t| write::mol_to_smiles(&t.mol, &RdkitSmilesParams::default()))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(RdkitTautomerEnumeration {
        smiles,
        keys: res.tautomers.keys().cloned().collect(),
        status: res.status,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tautomer_cancel_before_first_product() {
        let mol = crate::parse("CC=O").unwrap();
        let mut trace = Vec::new();
        let result = rdkit_enumerate_tautomers_with_callback(&mol, 1000, 1000, &mut |p| {
            trace.push(p.clone());
            false
        })
        .unwrap();
        assert_eq!(result.status, RdkitTautomerStatus::Canceled);
        assert_eq!(result.smiles, ["CC=O"]);
        assert_eq!(
            trace,
            [RdkitTautomerProgress {
                transforms: 1,
                keys: vec!["CC=O".into()]
            }]
        );
    }

    #[test]
    fn tautomer_limits_precede_callback() {
        let mol = crate::parse("CC=O").unwrap();
        for (max_tautomers, max_transforms, status) in [
            (1, 1000, RdkitTautomerStatus::MaxTautomersReached),
            (1000, 1, RdkitTautomerStatus::MaxTransformsReached),
            (1, 1, RdkitTautomerStatus::MaxTransformsReached),
        ] {
            let result = rdkit_enumerate_tautomers_with_callback(
                &mol,
                max_tautomers,
                max_transforms,
                &mut |_| panic!("limit must win"),
            )
            .unwrap();
            assert_eq!(result.status, status);
        }
        assert!(rdkit_enumerate_tautomers_with_callback(&mol, 0, 1, &mut |_| true).is_err());
        assert!(rdkit_enumerate_tautomers_with_callback(&mol, 1, 0, &mut |_| true).is_err());
    }

    #[test]
    fn tautomer_continue_equals_existing_path() {
        for input in ["CC=O", "O=c1cccc[nH]1", "NC(=N)N", "CC(=O)CC=O", "CCO", ""] {
            let mol = if input.is_empty() {
                chematic_core::MoleculeBuilder::new().build()
            } else {
                crate::parse(input).unwrap()
            };
            let expected = rdkit_enumerate_tautomers(&mol).unwrap();
            let result =
                rdkit_enumerate_tautomers_with_callback(&mol, 1000, 1000, &mut |_| true).unwrap();
            assert_eq!(result, expected, "{input}");
        }
    }
}
