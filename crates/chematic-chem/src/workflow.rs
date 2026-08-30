//! High-level, MCP-ready workflows built from the lower-level chematic APIs.
//!
//! These functions intentionally keep I/O out of the core logic. They return
//! structured Rust values that can be serialized by a future MCP/CLI/WASM layer.

use std::fmt;

use chematic_core::{Atom, AtomIdx, BondOrder, Element, Molecule, MoleculeBuilder};
use chematic_fp::{atom_pair_fp, ecfp4, maccs};
use chematic_smarts::{
    AtomPrimitive, AtomQuery, BondPrimitive, BondQuery, McsConfig, QueryMolecule,
    find_mcs_with_config,
};

use crate::{
    bertz_ct, detect_named_functional_groups, exact_mass, formal_charge_sum, fsp3, ghose_passes,
    hbd_count, heavy_atom_count, identify_functional_groups, labute_asa, logp_and_mr,
    molecular_weight, murcko_scaffold, num_heteroatoms, num_stereocenters,
    pains_passes_and_matches, qed_with_bundle, reos_passes, ring_bundle, sa_score_with_bundle,
    tpsa, wiener_index,
};

/// Error type for high-level workflow APIs.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum WorkflowError {
    /// A SMILES input failed to parse.
    SmilesParse {
        /// Index of the input in the caller-provided list.
        index: usize,
        /// Original SMILES string.
        smiles: String,
        /// Parser error message.
        message: String,
    },
    /// A molecule exceeded the configured atom-count limit.
    TooManyAtoms {
        /// Index of the input in the caller-provided list.
        index: usize,
        /// Original SMILES string.
        smiles: String,
        /// Actual atom count.
        atom_count: usize,
        /// Configured limit.
        max_atoms: usize,
    },
    /// The caller provided too many molecules for the requested workflow.
    TooManyMolecules {
        /// Actual input count.
        count: usize,
        /// Configured limit.
        max_molecules: usize,
    },
    /// The workflow requires at least two valid molecules.
    NeedAtLeastTwoMolecules,
}

impl fmt::Display for WorkflowError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::SmilesParse {
                index,
                smiles,
                message,
            } => {
                write!(
                    f,
                    "SMILES parse failed at index {index} ({smiles:?}): {message}"
                )
            }
            Self::TooManyAtoms {
                index,
                smiles,
                atom_count,
                max_atoms,
            } => {
                write!(
                    f,
                    "molecule at index {index} ({smiles:?}) has {atom_count} atoms; max is {max_atoms}"
                )
            }
            Self::TooManyMolecules {
                count,
                max_molecules,
            } => {
                write!(f, "input has {count} molecules; max is {max_molecules}")
            }
            Self::NeedAtLeastTwoMolecules => write!(f, "at least two molecules are required"),
        }
    }
}

impl std::error::Error for WorkflowError {}

/// Safety limits for LLM/tool-facing workflows.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct WorkflowLimits {
    /// Maximum molecules accepted by multi-molecule workflows.
    pub max_molecules: usize,
    /// Maximum atom count accepted for any single molecule.
    pub max_atoms: usize,
}

impl Default for WorkflowLimits {
    fn default() -> Self {
        Self {
            max_molecules: 256,
            max_atoms: 512,
        }
    }
}

/// Scalar descriptors and structural counts for a molecule.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct DescriptorSummary {
    pub molecular_weight: f64,
    pub exact_mass: f64,
    pub tpsa: f64,
    pub logp: f64,
    pub molar_refractivity: f64,
    pub hbd: usize,
    pub hba: usize,
    pub rotatable_bonds: usize,
    pub heavy_atom_count: usize,
    pub ring_count: usize,
    pub num_heteroatoms: usize,
    pub num_stereocenters: usize,
    pub num_spiro_atoms: usize,
    pub num_bridgehead_atoms: usize,
    pub fsp3: f64,
    pub qed: f64,
    pub sa_score: f64,
    pub formal_charge_sum: i32,
    pub labute_asa: f64,
    pub bertz_ct: f64,
    pub wiener_index: f64,
}

/// Common oral drug-likeness and structural-alert filter outcomes.
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct FilterSummary {
    pub lipinski_passes: bool,
    pub veber_passes: bool,
    pub egan_passes: bool,
    pub ghose_passes: bool,
    pub reos_passes: bool,
    pub pains_passes: bool,
    pub pains_alerts: Vec<String>,
}

/// Functional-group label with atom indices.
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct FunctionalGroupSummary {
    pub name: String,
    pub atom_indices: Vec<usize>,
}

/// Named group label with atom indices.
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct NamedGroupSummary {
    pub name: String,
    pub atom_indices: Vec<usize>,
}

/// Complete single-molecule report for tool-facing integrations.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct MoleculeReport {
    pub input_smiles: String,
    pub canonical_smiles: String,
    pub formula: String,
    pub murcko_scaffold_smiles: Option<String>,
    pub descriptors: DescriptorSummary,
    pub filters: FilterSummary,
    pub functional_groups: Vec<FunctionalGroupSummary>,
    pub named_groups: Vec<NamedGroupSummary>,
}

/// Options for single-molecule reports.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct ReportOptions {
    pub limits: WorkflowLimits,
}

/// Similarity values for a pair of molecules.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SimilaritySummary {
    pub ecfp4_tanimoto: f64,
    pub maccs_tanimoto: f64,
    pub atom_pair_tanimoto: f64,
}

/// Pairwise comparison entry.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct PairwiseComparison {
    pub left_index: usize,
    pub right_index: usize,
    pub similarities: SimilaritySummary,
}

/// Descriptor delta between two molecules (`right - left`).
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct DescriptorDelta {
    pub left_index: usize,
    pub right_index: usize,
    pub molecular_weight: f64,
    pub exact_mass: f64,
    pub tpsa: f64,
    pub logp: f64,
    pub hbd: isize,
    pub hba: isize,
    pub rotatable_bonds: isize,
    pub qed: f64,
    pub sa_score: f64,
}

/// Multi-molecule comparison report.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct MoleculeComparison {
    pub reports: Vec<MoleculeReport>,
    pub pairwise: Vec<PairwiseComparison>,
    pub descriptor_deltas: Vec<DescriptorDelta>,
    pub mcs_smiles: Option<String>,
}

/// Options for multi-molecule comparison.
#[derive(Debug, Clone)]
pub struct CompareOptions {
    pub limits: WorkflowLimits,
    pub mcs_config: McsConfig,
}

impl Default for CompareOptions {
    fn default() -> Self {
        Self {
            limits: WorkflowLimits::default(),
            mcs_config: McsConfig {
                timeout_ms: Some(250),
                ..McsConfig::default()
            },
        }
    }
}

/// Per-record result for batch screening.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ScreeningRecord {
    pub input_index: usize,
    pub input_smiles: String,
    pub report: Option<MoleculeReport>,
    pub error: Option<String>,
}

/// Batch screening report.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ScreeningReport {
    pub records: Vec<ScreeningRecord>,
    /// Original input indices selected by MaxMin diversity picking.
    pub maxmin_picks: Vec<usize>,
    /// Original input indices grouped by Butina clustering.
    pub butina_clusters: Vec<Vec<usize>>,
}

/// Options for batch screening.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ScreenOptions {
    pub limits: WorkflowLimits,
    pub maxmin_pick_count: usize,
    pub butina_cutoff: f64,
}

impl Default for ScreenOptions {
    fn default() -> Self {
        Self {
            limits: WorkflowLimits::default(),
            maxmin_pick_count: 0,
            butina_cutoff: 0.65,
        }
    }
}

/// Build a complete report from a SMILES string using default options.
pub fn molecule_report(smiles: &str) -> Result<MoleculeReport, WorkflowError> {
    molecule_report_with_options(smiles, &ReportOptions::default())
}

/// Build a complete report from a SMILES string.
pub fn molecule_report_with_options(
    smiles: &str,
    options: &ReportOptions,
) -> Result<MoleculeReport, WorkflowError> {
    let mol = parse_checked(0, smiles, options.limits.max_atoms)?;
    Ok(report_for_molecule(smiles, &mol))
}

/// Compare two or more SMILES strings using default options.
pub fn compare_molecules(smiles: &[&str]) -> Result<MoleculeComparison, WorkflowError> {
    compare_molecules_with_options(smiles, &CompareOptions::default())
}

/// Compare two or more SMILES strings.
pub fn compare_molecules_with_options(
    smiles: &[&str],
    options: &CompareOptions,
) -> Result<MoleculeComparison, WorkflowError> {
    if smiles.len() < 2 {
        return Err(WorkflowError::NeedAtLeastTwoMolecules);
    }
    if smiles.len() > options.limits.max_molecules {
        return Err(WorkflowError::TooManyMolecules {
            count: smiles.len(),
            max_molecules: options.limits.max_molecules,
        });
    }

    let mols = parse_many_checked(smiles, options.limits.max_atoms)?;
    let reports = smiles
        .iter()
        .zip(mols.iter())
        .map(|(s, mol)| report_for_molecule(s, mol))
        .collect::<Vec<_>>();

    let mut pairwise = Vec::new();
    let mut descriptor_deltas = Vec::new();
    for i in 0..mols.len() {
        for j in (i + 1)..mols.len() {
            pairwise.push(PairwiseComparison {
                left_index: i,
                right_index: j,
                similarities: similarity_summary(&mols[i], &mols[j]),
            });
            descriptor_deltas.push(descriptor_delta(i, j, &reports[i], &reports[j]));
        }
    }

    let mol_refs = mols.iter().collect::<Vec<_>>();
    let qmol = find_mcs_with_config(&mol_refs, &options.mcs_config);
    let mcs_smiles = query_molecule_to_smiles(&qmol);

    Ok(MoleculeComparison {
        reports,
        pairwise,
        descriptor_deltas,
        mcs_smiles,
    })
}

/// Screen a SMILES list using default options. Invalid records are retained
/// as per-record errors instead of failing the whole batch.
pub fn screen_smiles(smiles: &[&str]) -> ScreeningReport {
    screen_smiles_with_options(smiles, &ScreenOptions::default())
}

/// Screen a SMILES list. Invalid records are retained as per-record errors
/// instead of failing the whole batch.
pub fn screen_smiles_with_options(smiles: &[&str], options: &ScreenOptions) -> ScreeningReport {
    let mut records = Vec::with_capacity(smiles.len());
    let mut valid_mols = Vec::new();
    let mut valid_original_indices = Vec::new();

    for (idx, smi) in smiles.iter().enumerate() {
        if idx >= options.limits.max_molecules {
            records.push(ScreeningRecord {
                input_index: idx,
                input_smiles: (*smi).to_string(),
                report: None,
                error: Some(format!(
                    "input index exceeds max_molecules={}",
                    options.limits.max_molecules
                )),
            });
            continue;
        }

        match parse_checked(idx, smi, options.limits.max_atoms) {
            Ok(mol) => {
                let report = report_for_molecule(smi, &mol);
                valid_original_indices.push(idx);
                valid_mols.push(mol);
                records.push(ScreeningRecord {
                    input_index: idx,
                    input_smiles: (*smi).to_string(),
                    report: Some(report),
                    error: None,
                });
            }
            Err(err) => records.push(ScreeningRecord {
                input_index: idx,
                input_smiles: (*smi).to_string(),
                report: None,
                error: Some(err.to_string()),
            }),
        }
    }

    let maxmin_picks = if options.maxmin_pick_count == 0 || valid_mols.is_empty() {
        Vec::new()
    } else {
        crate::maxmin_picks(&valid_mols, options.maxmin_pick_count, |a, b| {
            ecfp4(a).tanimoto(&ecfp4(b))
        })
        .into_iter()
        .map(|valid_idx| valid_original_indices[valid_idx])
        .collect()
    };

    let butina_clusters = if valid_mols.is_empty() {
        Vec::new()
    } else {
        crate::butina_cluster(&valid_mols, options.butina_cutoff, |a, b| {
            ecfp4(a).tanimoto(&ecfp4(b))
        })
        .into_iter()
        .map(|cluster| {
            cluster
                .into_iter()
                .map(|valid_idx| valid_original_indices[valid_idx])
                .collect()
        })
        .collect()
    };

    ScreeningReport {
        records,
        maxmin_picks,
        butina_clusters,
    }
}

fn parse_many_checked(smiles: &[&str], max_atoms: usize) -> Result<Vec<Molecule>, WorkflowError> {
    smiles
        .iter()
        .enumerate()
        .map(|(idx, smi)| parse_checked(idx, smi, max_atoms))
        .collect()
}

fn parse_checked(index: usize, smiles: &str, max_atoms: usize) -> Result<Molecule, WorkflowError> {
    let mol = chematic_smiles::parse(smiles).map_err(|e| WorkflowError::SmilesParse {
        index,
        smiles: smiles.to_string(),
        message: e.to_string(),
    })?;
    if mol.atom_count() > max_atoms {
        return Err(WorkflowError::TooManyAtoms {
            index,
            smiles: smiles.to_string(),
            atom_count: mol.atom_count(),
            max_atoms,
        });
    }
    Ok(mol)
}

fn report_for_molecule(input_smiles: &str, mol: &Molecule) -> MoleculeReport {
    let scaffold = murcko_scaffold(mol);
    let murcko_scaffold_smiles = if scaffold.atom_count() == 0 {
        None
    } else {
        Some(chematic_smiles::canonical_smiles(&scaffold))
    };

    // Compute all ring-derived descriptors with a single find_sssr call.
    let rb = ring_bundle(mol);
    let mw = molecular_weight(mol);
    // Compute LogP and MR together to share the 117-pattern Crippen SMARTS pass.
    let (logp, mr) = logp_and_mr(mol);
    let tpsa_val = tpsa(mol);
    // Single explicit-H + SSSR + 480-pattern pass for both PAINS flag and alert names.
    let (pains_ok, pains_alert_names) = pains_passes_and_matches(mol);

    MoleculeReport {
        input_smiles: input_smiles.to_string(),
        canonical_smiles: chematic_smiles::canonical_smiles(mol),
        formula: mol.total_formula(),
        murcko_scaffold_smiles,
        descriptors: DescriptorSummary {
            molecular_weight: mw,
            exact_mass: exact_mass(mol),
            tpsa: tpsa_val,
            logp,
            molar_refractivity: mr,
            hbd: hbd_count(mol),
            hba: rb.hba_count,
            rotatable_bonds: rb.rotatable_bond_count,
            heavy_atom_count: heavy_atom_count(mol),
            ring_count: rb.ring_count,
            num_heteroatoms: num_heteroatoms(mol),
            num_stereocenters: num_stereocenters(mol),
            num_spiro_atoms: rb.num_spiro_atoms,
            num_bridgehead_atoms: rb.num_bridgehead_atoms,
            fsp3: fsp3(mol),
            qed: qed_with_bundle(mol, &rb),
            sa_score: sa_score_with_bundle(mol, &rb),
            formal_charge_sum: formal_charge_sum(mol),
            labute_asa: labute_asa(mol),
            bertz_ct: bertz_ct(mol),
            wiener_index: wiener_index(mol),
        },
        filters: FilterSummary {
            // Inline only the filters that call rotatable_bond_count or hba_count (→ find_sssr).
            // Other filters (egan, ghose, reos) don't use ring data — call them directly.
            lipinski_passes: mw <= 500.0
                && hbd_count(mol) <= 5
                && rb.hba_count <= 10
                && logp <= 5.0,
            veber_passes: rb.rotatable_bond_count <= 10 && tpsa_val <= 140.0,
            egan_passes: tpsa_val <= 131.6 && logp <= 5.88,
            ghose_passes: ghose_passes(mol),
            reos_passes: reos_passes(mol),
            pains_passes: pains_ok,
            pains_alerts: pains_alert_names.into_iter().map(str::to_string).collect(),
        },
        functional_groups: identify_functional_groups(mol)
            .into_iter()
            .map(|fg| FunctionalGroupSummary {
                name: fg.atom_types,
                atom_indices: fg.atom_indices,
            })
            .collect(),
        named_groups: detect_named_functional_groups(mol)
            .into_iter()
            .map(|ng| NamedGroupSummary {
                name: ng.name.to_string(),
                atom_indices: ng.atoms.into_iter().map(|idx| idx.0 as usize).collect(),
            })
            .collect(),
    }
}

fn similarity_summary(a: &Molecule, b: &Molecule) -> SimilaritySummary {
    SimilaritySummary {
        ecfp4_tanimoto: ecfp4(a).tanimoto(&ecfp4(b)),
        maccs_tanimoto: maccs(a).tanimoto(&maccs(b)),
        atom_pair_tanimoto: atom_pair_fp(a).tanimoto(&atom_pair_fp(b)),
    }
}

fn descriptor_delta(
    left_index: usize,
    right_index: usize,
    left: &MoleculeReport,
    right: &MoleculeReport,
) -> DescriptorDelta {
    DescriptorDelta {
        left_index,
        right_index,
        molecular_weight: right.descriptors.molecular_weight - left.descriptors.molecular_weight,
        exact_mass: right.descriptors.exact_mass - left.descriptors.exact_mass,
        tpsa: right.descriptors.tpsa - left.descriptors.tpsa,
        logp: right.descriptors.logp - left.descriptors.logp,
        hbd: right.descriptors.hbd as isize - left.descriptors.hbd as isize,
        hba: right.descriptors.hba as isize - left.descriptors.hba as isize,
        rotatable_bonds: right.descriptors.rotatable_bonds as isize
            - left.descriptors.rotatable_bonds as isize,
        qed: right.descriptors.qed - left.descriptors.qed,
        sa_score: right.descriptors.sa_score - left.descriptors.sa_score,
    }
}

fn query_molecule_to_smiles(qmol: &QueryMolecule) -> Option<String> {
    if qmol.atoms.is_empty() {
        return None;
    }

    // `build_query`/`molecule_to_query` (chematic-smarts) always wrap each
    // atom's query as `And(AtomicNum(n), Aromatic(bool))`, never a bare
    // `AtomicNum` primitive, so this must recurse through `And` or it
    // silently falls through to Carbon for every atom.
    fn extract_atomic_num(q: &AtomQuery) -> Option<u8> {
        match q {
            AtomQuery::Primitive(AtomPrimitive::AtomicNum(n)) => Some(*n),
            AtomQuery::And(lhs, rhs) => extract_atomic_num(lhs).or_else(|| extract_atomic_num(rhs)),
            _ => None,
        }
    }

    let mut aromatic_atoms = vec![false; qmol.atoms.len()];
    for (atom_idx, neighbors) in qmol.adj.iter().enumerate() {
        for (bond_idx, neighbor_idx) in neighbors {
            if matches!(
                qmol.bonds[*bond_idx].query,
                BondQuery::Primitive(BondPrimitive::Aromatic)
            ) {
                aromatic_atoms[atom_idx] = true;
                aromatic_atoms[*neighbor_idx] = true;
            }
        }
    }

    let mut builder = MoleculeBuilder::new();
    for (idx, qa) in qmol.atoms.iter().enumerate() {
        let elem = extract_atomic_num(&qa.query)
            .and_then(Element::from_atomic_number)
            .unwrap_or(Element::C);
        let mut atom = Atom::new(elem);
        atom.aromatic = aromatic_atoms[idx];
        builder.add_atom(atom);
    }

    for (atom_idx, neighbors) in qmol.adj.iter().enumerate() {
        for (bond_idx, neighbor_idx) in neighbors {
            if atom_idx < *neighbor_idx {
                let order = match &qmol.bonds[*bond_idx].query {
                    BondQuery::Primitive(BondPrimitive::Double) => BondOrder::Double,
                    BondQuery::Primitive(BondPrimitive::Triple) => BondOrder::Triple,
                    BondQuery::Primitive(BondPrimitive::Aromatic) => BondOrder::Aromatic,
                    _ => BondOrder::Single,
                };
                let _ = builder.add_bond(
                    AtomIdx(atom_idx as u32),
                    AtomIdx(*neighbor_idx as u32),
                    order,
                );
            }
        }
    }

    Some(chematic_smiles::canonical_smiles(&builder.build()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn molecule_report_aspirin_has_core_fields() {
        let report = molecule_report("CC(=O)Oc1ccccc1C(=O)O").unwrap();

        assert_eq!(report.formula, "C9H8O4");
        assert_eq!(report.descriptors.heavy_atom_count, 13);
        assert!(report.descriptors.molecular_weight > 180.0);
        assert!(report.descriptors.tpsa > 60.0);
        assert!(report.filters.lipinski_passes);
        assert_eq!(report.murcko_scaffold_smiles.as_deref(), Some("c1ccccc1"));
    }

    #[test]
    fn molecule_report_invalid_smiles_returns_structured_error() {
        let err = molecule_report("C1CC").unwrap_err();
        assert!(matches!(err, WorkflowError::SmilesParse { index: 0, .. }));
    }

    #[test]
    fn compare_molecules_returns_pairwise_and_mcs() {
        let comparison = compare_molecules(&["c1ccccc1", "Cc1ccccc1"]).unwrap();

        assert_eq!(comparison.reports.len(), 2);
        assert_eq!(comparison.pairwise.len(), 1);
        assert_eq!(comparison.descriptor_deltas.len(), 1);
        assert!(comparison.pairwise[0].similarities.ecfp4_tanimoto > 0.0);
        assert_eq!(comparison.mcs_smiles.as_deref(), Some("c1ccccc1"));
    }

    #[test]
    fn compare_molecules_mcs_preserves_heteroatoms_and_aromaticity() {
        // Regression test: `query_molecule_to_smiles` once matched only the bare
        // `AtomQuery::Primitive(AtomicNum(n))` case, but at the time `find_mcs`'s
        // query atoms were always the compound `And(AtomicNum(n), Aromatic(bool))`
        // -- so every atom silently fell through to carbon. Benzene/toluene's
        // all-carbon MCS (see the test above) could never catch this; a
        // hydroxyphenyl MCS shared between aspirin and paracetamol can, since it
        // contains an O atom.
        let comparison =
            compare_molecules(&["CC(=O)Oc1ccccc1C(=O)O", "CC(=O)Nc1ccc(O)cc1"]).unwrap();
        let mcs = comparison
            .mcs_smiles
            .expect("aspirin/paracetamol share a hydroxyphenyl MCS");
        assert!(
            mcs.contains('O') || mcs.contains('o'),
            "MCS lost its oxygen atom: {mcs:?}"
        );
    }

    #[test]
    fn compare_molecules_needs_two_inputs() {
        let err = compare_molecules(&["CC"]).unwrap_err();
        assert_eq!(err, WorkflowError::NeedAtLeastTwoMolecules);
    }

    #[test]
    fn screen_smiles_keeps_invalid_records_and_original_indices() {
        let options = ScreenOptions {
            maxmin_pick_count: 2,
            ..ScreenOptions::default()
        };
        let report = screen_smiles_with_options(&["CC", "C1CC", "c1ccccc1"], &options);

        assert_eq!(report.records.len(), 3);
        assert!(report.records[0].report.is_some());
        assert!(report.records[1].report.is_none());
        assert!(report.records[1].error.is_some());
        assert!(report.records[2].report.is_some());
        assert!(report.maxmin_picks.iter().all(|idx| *idx == 0 || *idx == 2));
        assert!(
            report
                .butina_clusters
                .iter()
                .flatten()
                .all(|idx| *idx == 0 || *idx == 2)
        );
    }

    #[test]
    fn limits_reject_large_molecule() {
        let options = ReportOptions {
            limits: WorkflowLimits {
                max_atoms: 1,
                ..WorkflowLimits::default()
            },
        };
        let err = molecule_report_with_options("CC", &options).unwrap_err();
        assert!(matches!(
            err,
            WorkflowError::TooManyAtoms { atom_count: 2, .. }
        ));
    }

    #[test]
    fn molecule_report_sulfur_compound() {
        let report = molecule_report("c1ccccc1S(=O)(=O)N").unwrap();
        assert_eq!(report.descriptors.num_heteroatoms, 4); // S, O, O, N
        assert!(report.descriptors.molecular_weight > 150.0);
    }

    #[test]
    fn molecule_report_halogenated() {
        let report = molecule_report("ClC(Br)(F)I").unwrap();
        assert_eq!(report.descriptors.heavy_atom_count, 5); // C, Cl, Br, F, I
        assert!(report.descriptors.num_heteroatoms == 4); // Cl, Br, F, I (not C)
    }

    #[test]
    fn molecule_report_complex_aromatic() {
        // Quinoline (bicyclic aromatic)
        let report = molecule_report("c1ccc2ncccc2c1").unwrap();
        assert_eq!(report.descriptors.heavy_atom_count, 10); // 9 C + 1 N
        assert!(report.descriptors.ring_count >= 2);
        assert!(report.filters.lipinski_passes);
    }

    #[test]
    fn molecule_report_large_valid_molecule() {
        // Taxol-like structure (simplified, ~50 atoms)
        let report = molecule_report("CC(=O)Oc1ccccc1C(=O)N[C@@H]1C[C@H]2CC(C)(C)[C@@H](O)C[C@]2(OC(=O)C(C)C)[C@]1(O)C(=O)c1ccccc1").unwrap();
        assert!(report.descriptors.heavy_atom_count > 30);
        assert!(report.descriptors.num_stereocenters > 0);
    }

    #[test]
    fn molecule_report_charge_species() {
        let report = molecule_report("[Na+].[Cl-]").unwrap();
        assert_eq!(report.formula, "ClNa");
        assert_eq!(report.descriptors.formal_charge_sum, 0); // +1 + -1 = 0
    }

    #[test]
    fn compare_molecules_identical_molecules() {
        let comparison = compare_molecules(&["c1ccccc1", "c1ccccc1"]).unwrap();
        let sim = comparison.pairwise[0].similarities.ecfp4_tanimoto;
        assert!(
            (sim - 1.0).abs() < 1e-6,
            "identical molecules should have ~100% similarity"
        );
    }

    #[test]
    fn screen_smiles_empty_batch() {
        let report = screen_smiles(&[]);
        assert_eq!(report.records.len(), 0);
        assert_eq!(report.maxmin_picks.len(), 0);
    }

    #[test]
    fn molecule_report_aromatic_nitrogen() {
        // Pyridine
        let report = molecule_report("c1ccncc1").unwrap();
        assert_eq!(report.descriptors.ring_count, 1);
        assert!(report.descriptors.hba > 0); // N is an acceptor
    }
}
