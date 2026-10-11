// ---------------------------------------------------------------------------
// 8. LogP (Wildman-Crippen, calibrated)
// ---------------------------------------------------------------------------

/// Crippen types on the RDKit port's AddHs molecule (diagnostics).
#[doc(hidden)]
pub fn crippen_types_rdkit_model(mol: &Molecule) -> Option<Vec<Option<(String, f64)>>> {
    let patterns: Vec<&str> = CRIPPEN_SMARTS.iter().map(|e| e.0).collect();
    let (types, _) = chematic_smiles::rdkit_addhs_first_pattern(mol, &patterns)?;
    Some(
        types
            .into_iter()
            .map(|t| t.map(|t| (CRIPPEN_SMARTS[t].0.to_string(), CRIPPEN_SMARTS[t].1)))
            .collect(),
    )
}

/// Whether Crippen typing must run on the RDKit port's molecule: chematic's
/// perception may differ from RDKit's (see
/// [`chematic_perception::rdkit_model_may_disagree`]), hydrogen graph atoms,
/// or charged hydroxyl/oxonium hydrogens (`[OH-]`, `[OH+]`).
fn crippen_needs_rdkit_model(mol: &Molecule) -> bool {
    chematic_perception::rdkit_model_may_disagree(mol)
        || mol.atoms().any(|(_, a)| {
            !a.wildcard
                && (a.element == Element::H
                    || (a.element == Element::O
                        && a.charge != 0
                        && a.hydrogen_count.is_some_and(|h| h > 0)))
        })
}

/// [`crippen_totals`] typed on the RDKit port's `AddHs(MolFromSmiles(s))`.
fn crippen_totals_rdkit_model(mol: &Molecule) -> Option<(f64, f64)> {
    chematic_smiles::rdkit_crippen_logp_mr(mol)
}

/// Compute a Wildman-Crippen LogP (RDKit-compatible).
///
/// Uses the same 117-entry SMARTS atom-type table as RDKit `Crippen.MolLogP()`,
/// taken directly from `rdkit/Chem/Crippen.py` (Wildman & Crippen 1999). The
/// first matching pattern for each atom determines its contribution. Results are
/// RDKit-compatible for common drug-like molecules.
///
/// RDKit Crippen LogP atom-type table (Wildman & Crippen 1999).
///
/// Ordered exactly as in RDKit `Crippen.py`. First matching pattern for each atom wins.
/// Each entry: (smarts, logp_contribution, mr_contribution).
/// Patterns that fail to parse are silently skipped.
static CRIPPEN_SMARTS: &[(&str, f64, f64)] = &[
    // --- Carbon (C1-C27 + CS fallback) ---
    ("[CH4]", 0.1441, 2.503),
    ("[CH3]C", 0.1441, 2.503),
    ("[CH2](C)C", 0.1441, 2.503),
    ("[CH](C)(C)C", 0.0000, 2.433),
    ("[C](C)(C)(C)C", 0.0000, 2.433),
    ("[CH3][N,O,P,S,F,Cl,Br,I]", -0.2035, 2.753),
    ("[CH2X4]([N,O,P,S,F,Cl,Br,I])[A;!#1]", -0.2035, 2.753),
    (
        "[CH1X4]([N,O,P,S,F,Cl,Br,I])([A;!#1])[A;!#1]",
        -0.2051,
        2.731,
    ),
    (
        "[CH0X4]([N,O,P,S,F,Cl,Br,I])([A;!#1])([A;!#1])[A;!#1]",
        -0.2051,
        2.731,
    ),
    ("[C]=[!C;A;!#1]", -0.2783, 5.007),
    ("[CH2]=C", 0.1551, 3.513),
    ("[CH1](=C)[A;!#1]", 0.1551, 3.513),
    ("[CH0](=C)([A;!#1])[A;!#1]", 0.1551, 3.513),
    ("[C](=C)=C", 0.1551, 3.513),
    ("[CX2]#[A;!#1]", 0.0017, 3.888),
    ("[CH3]c", 0.08452, 2.464),
    ("[CH3]a", -0.1444, 2.412),
    ("[CH2X4]a", -0.0516, 2.488),
    ("[CHX4]a", 0.1193, 2.582),
    ("[CH0X4]a", -0.0967, 2.576),
    ("[cH0]-[A;!C;!N;!O;!S;!F;!Cl;!Br;!I;!#1]", -0.5443, 4.041),
    ("[c][#9]", 0.0000, 3.257),
    ("[c][#17]", 0.2450, 3.564),
    ("[c][#35]", 0.1980, 3.180),
    ("[c][#53]", 0.0000, 3.104),
    ("[cH]", 0.1581, 3.350),
    ("[c](:a)(:a):a", 0.2955, 4.346),
    ("[c](:a)(:a)-a", 0.2713, 3.904),
    ("[c](:a)(:a)-C", 0.1360, 3.509),
    ("[c](:a)(:a)-N", 0.4619, 4.067),
    ("[c](:a)(:a)-O", 0.5437, 3.853),
    ("[c](:a)(:a)-S", 0.1893, 2.673),
    ("[c](:a)(:a)=[C,N,O]", -0.8186, 3.135),
    ("[C](=C)(a)[A;!#1]", 0.2640, 4.305),
    ("[C](=C)(c)a", 0.2640, 4.305),
    ("[CH1](=C)a", 0.2640, 4.305),
    ("[C]=c", 0.2640, 4.305),
    ("[CX4][A;!C;!N;!O;!P;!S;!F;!Cl;!Br;!I;!#1]", 0.2148, 2.693),
    ("[#6]", 0.08129, 3.243), // CS fallback
    // --- Hydrogen (H1-H4 + HS fallback) ---
    ("[#1][#6,#1]", 0.1230, 1.057),
    ("[#1]O[CX4,c]", -0.2677, 1.395),
    ("[#1]O[!C;!N;!O;!S]", -0.2677, 1.395),
    ("[#1][!C;!N;!O]", -0.2677, 1.395),
    ("[#1][#7]", 0.2142, 0.9627),
    ("[#1]O[#7]", 0.2142, 0.9627),
    ("[#1]OC=[#6,#7,O,S]", 0.2980, 1.805),
    ("[#1]O[O,S]", 0.2980, 1.805),
    ("[#1]", 0.1125, 1.112), // HS fallback
    // --- Nitrogen (N1-N14 + NS fallback) ---
    ("[NH2+0][A;!#1]", -1.0190, 2.262),
    ("[NH+0]([A;!#1])[A;!#1]", -0.7096, 2.173),
    ("[NH2+0]a", -1.0270, 2.827),
    ("[NH1+0]([!#1;A,a])a", -0.5188, 3.000),
    ("[NH+0]=[!#1;A,a]", 0.08387, 1.757),
    ("[N+0](=[!#1;A,a])[!#1;A,a]", 0.1836, 2.428),
    ("[N+0]([A;!#1])([A;!#1])[A;!#1]", -0.3187, 1.839),
    ("[N+0](a)([!#1;A,a])[A;!#1]", -0.4458, 2.819),
    ("[N+0](a)(a)a", -0.4458, 2.819),
    ("[N+0]#[A;!#1]", 0.01508, 1.725),
    ("[NH3,NH2,NH;+,+2,+3]", -1.9500, 0.000),
    ("[n+0]", -0.3239, 2.202),
    ("[n;+,+2,+3]", -1.1190, 0.000),
    (
        "[NH0;+,+2,+3]([A;!#1])([A;!#1])([A;!#1])[A;!#1]",
        -0.3396,
        0.2604,
    ),
    ("[NH0;+,+2,+3](=[A;!#1])([A;!#1])[!#1;A,a]", -0.3396, 0.2604),
    ("[NH0;+,+2,+3](=[#6])=[#7]", -0.3396, 0.2604),
    ("[N;+,+2,+3]#[A;!#1]", 0.2887, 3.359),
    ("[N;-,-2,-3]", 0.2887, 3.359),
    ("[N;+,+2,+3](=[N;-,-2,-3])=N", 0.2887, 3.359),
    ("[#7]", -0.4806, 2.134), // NS fallback
    // --- Oxygen (O1-O12 + OS fallback) ---
    ("[o]", 0.1552, 1.080),
    ("[OH,OH2]", -0.2893, 0.8238),
    ("[O]([A;!#1])[A;!#1]", -0.0684, 1.085),
    ("[O](a)[!#1;A,a]", -0.4195, 1.182),
    ("[O]=[#7,#8]", 0.0335, 3.367),
    ("[OX1;-,-2,-3][#7]", 0.0335, 3.367),
    ("[OX1;-,-2,-3][#16]", -0.3339, 0.7774),
    ("[O;-0]=[#16;-0]", -0.3339, 0.7774),
    ("[O-]C(=O)", -1.3260, 0.000),
    ("[OX1;-,-2,-3][!#1;!N;!S]", -1.1890, 0.000),
    ("[O]=c", 0.1788, 3.135),
    ("[O]=[CH]C", -0.1526, 0.000),
    ("[O]=C(C)([A;!#1])", -0.1526, 0.000),
    ("[O]=[CH][N,O]", -0.1526, 0.000),
    ("[O]=[CH2]", -0.1526, 0.000),
    ("[O]=[CX2]=O", -0.1526, 0.000),
    ("[O]=[CH]c", 0.1129, 0.2215),
    ("[O]=C([C,c])[a;!#1]", 0.1129, 0.2215),
    ("[O]=C(c)[A;!#1]", 0.1129, 0.2215),
    ("[O]=C([!#1;!#6])[!#1;!#6]", 0.4833, 0.3890),
    ("[#8]", -0.1188, 0.6865), // OS fallback
    // --- Sulfur (S1-S3) ---
    ("[S;-,-2,-3,-4,+1,+2,+3,+5,+6]", -0.0024, 7.365),
    ("[S;-0]=[N,O,P,S]", -0.0024, 7.365),
    ("[S;A]", 0.6482, 7.591),
    ("[s;a]", 0.6237, 6.691),
    // --- Phosphorus ---
    ("[#15]", 0.8612, 6.920),
    // --- Halogens ---
    ("[#9;-0]", 0.4202, 1.108),
    ("[#17;-0]", 0.6895, 5.853),
    ("[#35;-0]", 0.8456, 8.927),
    ("[#53;-0]", 0.8857, 14.020),
    ("[#9,#17,#35,#53;-]", -2.9960, 0.000),
    ("[#53;+,+2,+3]", -2.9960, 0.000),
    ("[+;#3,#11,#19,#37,#55]", -2.9960, 0.000),
    // --- Metals ---
    ("[#3,#11,#19,#37,#55]", -0.3808, 5.754),
    ("[#4,#12,#20,#38,#56]", -0.3808, 5.754),
    ("[#5,#13,#31,#49,#81]", -0.3808, 5.754),
    ("[#14,#32,#50,#82]", -0.3808, 5.754),
    ("[#33,#51,#83]", -0.3808, 5.754),
    ("[#34,#52,#84]", -0.3808, 5.754),
    ("[#21,#22,#23,#24,#25,#26,#27,#28,#29,#30]", -0.0025, 0.000),
    ("[#39,#40,#41,#42,#43,#44,#45,#46,#47,#48]", -0.0025, 0.000),
    ("[#72,#73,#74,#75,#76,#77,#78,#79,#80]", -0.0025, 0.000),
];

type CrippenQueries = Vec<(Option<chematic_smarts::QueryMolecule>, f64, f64)>;
static CRIPPEN_QUERIES: OnceLock<CrippenQueries> = OnceLock::new();

fn get_crippen_queries() -> &'static CrippenQueries {
    CRIPPEN_QUERIES.get_or_init(|| {
        CRIPPEN_SMARTS
            .iter()
            .map(|&(sma, lp, mr)| (parse_smarts(sma).ok(), lp, mr))
            .collect()
    })
}

/// Pre-compute for each SMARTS pattern which target atoms satisfy query-atom-0.
/// Amortizes VF2 cost from O(n_atoms × n_patterns) → O(n_patterns) per molecule.
/// SSSR is computed once and shared across all 117 Crippen patterns.
fn crippen_anchor_types(
    mol: &Molecule,
    queries: &CrippenQueries,
) -> std::sync::Arc<Vec<Option<usize>>> {
    // For each atom, the first Crippen pattern (table order) with an embedding
    // anchored at that atom as query atom 0. Equivalent to the historical
    // "collect every non-uniquified match of every pattern, then take the
    // first pattern whose anchor set contains the atom", but stops at the
    // first embedding. Memoized on the (already memoized) aromatic view, so
    // LogP and MR share one typing pass.
    mol.derived(chematic_core::DerivedSlot::CrippenTypes, || {
        let rings = chematic_perception::find_sssr_shared(mol);
        let refs: Vec<Option<&chematic_smarts::QueryMolecule>> =
            queries.iter().map(|(q, _, _)| q.as_ref()).collect();
        chematic_smarts::first_anchored_match_per_atom(&refs, mol, &rings)
    })
}

/// Wildman-Crippen LogP per-atom contributions (RDKit-compatible SMARTS dispatch).
///
/// Uses the same 117-entry atom-type table as RDKit `Crippen.MolLogP()`. Each atom
/// is assigned the contribution of the first matching SMARTS pattern. H contributions
/// are included via explicit `[#1]` patterns; heavy atoms and their implicit H are
/// handled together as in the original Wildman-Crippen definition.
/// Index matches mol.atoms().
pub fn logp_crippen_per_atom(mol: &Molecule) -> Vec<f64> {
    let queries = get_crippen_queries();
    // RDKit assigns Crippen atom types after aromaticity perception. Match on
    // that representation so an equivalent Kekulé spelling receives the same
    // aromatic SMARTS type; contributions stay indexed to the caller's graph.
    let mol_arom = descriptor_aromaticity(mol);
    // Pre-compute once per molecule: for each pattern, which atoms satisfy query-atom-0?
    // Previously O(n_atoms × n_patterns × VF2); now O(n_patterns × VF2 + n_atoms × n_patterns).
    let anchor_types = crippen_anchor_types(&mol_arom, queries);

    let h_fallback = CRIPPEN_SMARTS
        .iter()
        .find(|(sma, _, _)| *sma == "[#1]")
        .map(|e| e.1)
        .unwrap_or(0.1125);

    mol.atoms()
        .map(|(idx, atom)| {
            if atom.element.atomic_number() == 1 {
                return 0.0;
            }
            // Implicit H from valence rules + explicit H neighbors (e.g. [H], [2H], [3H]).
            // Explicit isotopic H in the graph contribute 0.0 themselves but their LogP
            // must be counted via their parent heavy atom (same as regular implicit H).
            let h_count_impl = implicit_hcount(mol, idx);
            let h_count_expl = mol
                .neighbors(idx)
                .filter(|(nb, _)| mol.atom(*nb).element.atomic_number() == 1)
                .count() as u8;
            let h_count = h_count_impl + h_count_expl;
            let heavy = anchor_types[idx.0 as usize]
                .map(|t| queries[t].1)
                .unwrap_or(0.0);
            // RDKit adds each hydrogen's contribution to its heavy atom one
            // at a time (getCrippenAtomContribs folds the added H atoms
            // back); a single multiplication rounds differently.
            let mut contrib = heavy;
            if h_count > 0 {
                let h = h_logp_for_parent(
                    &mol_arom,
                    idx,
                    atom.element.atomic_number(),
                    // Hydrogen typing sees the RDKit-perceived environment, like
                    // the heavy-atom SMARTS typing (Kekule phenol OH is H2, not H4).
                    mol_arom.atom(idx).aromatic,
                    h_fallback,
                );
                for _ in 0..h_count {
                    contrib += h;
                }
            }
            contrib
        })
        .collect()
}

/// Whether a hydrogen on the (non-aromatic) oxygen `parent` is Crippen type
/// H2 by one of RDKit's first-listed patterns, `[#1]O[CX4,c]` or
/// `[#1]O[!#6;!#7;!#8;!#16]`. RDKit takes the first matching type in table
/// order, so these win over H3 (`[#1]O[#7]`) and H4 (`[#1]OC=[#6,#7,O,S]`)
/// when the oxygen has several neighbours, e.g. a metal-bound `[OH+]`
/// of a carboxylate.
fn oh_hydrogen_is_h2(mol: &Molecule, parent: AtomIdx) -> bool {
    mol.neighbors(parent).any(|(nb, _)| {
        let a = mol.atom(nb);
        let z = if a.wildcard {
            0
        } else {
            a.element.atomic_number()
        };
        match z {
            6 => a.aromatic || mol.degree(nb) + usize::from(implicit_hcount(mol, nb)) == 4,
            7 | 8 | 16 => false,
            _ => true,
        }
    })
}

/// Determine the LogP contribution per implicit-H attached to `parent_idx`.
/// Matches against H-type SMARTS by looking at the parent atom's environment.
fn h_logp_for_parent(
    mol: &Molecule,
    parent_idx: AtomIdx,
    parent_an: u8,
    parent_aromatic: bool,
    fallback: f64,
) -> f64 {
    // H type depends on what atom the H is bonded to:
    // H1: [#1][#6,#1]         → H on C  → 0.1230
    // H2: [#1]O[CX4,c] etc.  → H on O (various)  → -0.2677
    // H3: [#1][#7]            → H on N  → 0.2142
    // H4: [#1]OC=[#6,#7,O,S] → H on O adjacent to C=O → 0.2980
    // HS: [#1]                → fallback → 0.1125
    //
    // Since we can't easily run SMARTS with a virtual H atom, we use the
    // original hand-coded H dispatch which is already accurate:
    match parent_an {
        6 => 0.1230, // H1: H on C
        7 => 0.2142, // H3: H on N
        8 => {
            if parent_aromatic {
                fallback // aromatic O (furan-type) — HS fallback
            } else if oh_hydrogen_is_h2(mol, parent_idx) {
                -0.2677
            } else if mol.neighbors(parent_idx).any(|(nb, _)| {
                // H4: O bonded to C, and that C has a double bond to C/N/O/S
                // Matches RDKit [#1]OC=[#6,#7,O,S]: covers carboxylic, enol, vinylogous OH
                mol.atom(nb).element.atomic_number() == 6
                    && mol.neighbors(nb).any(|(nb2, bidx2)| {
                        nb2 != parent_idx
                            && mol.bond(bidx2).order == BondOrder::Double
                            && matches!(mol.atom(nb2).element.atomic_number(), 6 | 7 | 8 | 16)
                    })
            }) {
                // H4: H on O adjacent to C=X (X = C, N, O, S)
                0.2980
            } else if mol.neighbors(parent_idx).any(|(nb, _)| {
                mol.atom(nb).element.atomic_number() == 7 // O bonded to N
            }) {
                // [#1]O[#7]: H on O attached to N — oxime / hydroxamic acid
                0.2142
            } else if mol.neighbors(parent_idx).any(|(nb, _)| {
                let an = mol.atom(nb).element.atomic_number();
                an == 8 || an == 16 // [#1]O[O,S]: peroxide / sulfonic / sulfuric acid OH
            }) {
                0.2980
            } else {
                // H2: [#1]O[CX4,c] — aliphatic and phenolic OH both use -0.2677
                -0.2677
            }
        }
        // H2: [#1][!C;!N;!O] — H on S, P, Si, halogens, etc.
        _ => -0.2677,
    }
}

/// RDKit's `getCrippenAtomContribs` on the molecule as given (no added
/// hydrogens): each atom's own type contribution `(logP, MR)`, hydrogens
/// not folded in. A graph hydrogen gets the hydrogen type of its neighbour.
/// The SlogP/SMR VSA descriptors bin atoms by these values.
pub(crate) fn crippen_atom_type_contribs(mol: &Molecule) -> Vec<(f64, f64)> {
    if crippen_needs_rdkit_model(mol)
        && let Some(c) = chematic_smiles::rdkit_crippen_contribs_no_hs(mol)
        && c.len() == mol.atom_count()
    {
        return c;
    }
    let queries = get_crippen_queries();
    let mol_arom = descriptor_aromaticity(mol);
    let anchor_types = crippen_anchor_types(&mol_arom, queries);
    let (h_logp_fallback, h_mr_fallback) = CRIPPEN_SMARTS
        .iter()
        .find(|(sma, _, _)| *sma == "[#1]")
        .map(|e| (e.1, e.2))
        .unwrap_or((0.1125, 1.112));
    mol.atoms()
        .map(|(idx, atom)| {
            if atom.element.atomic_number() == 1 && !atom.wildcard {
                return match mol
                    .neighbors(idx)
                    .find(|(nb, _)| mol.atom(*nb).element.atomic_number() != 1)
                {
                    Some((parent, _)) => {
                        let an = mol.atom(parent).element.atomic_number();
                        let aromatic = mol_arom.atom(parent).aromatic;
                        (
                            h_logp_for_parent(&mol_arom, parent, an, aromatic, h_logp_fallback),
                            h_mr_for_parent(&mol_arom, parent, an, aromatic, h_mr_fallback),
                        )
                    }
                    None => (h_logp_fallback, h_mr_fallback),
                };
            }
            anchor_types[idx.0 as usize]
                .map(|t| (queries[t].1, queries[t].2))
                .unwrap_or((0.0, 0.0))
        })
        .collect()
}

/// Compute the Crippen log P (octanol/water partition coefficient) of `mol`.
///
/// The per-atom contributions of [`logp_crippen_per_atom`], summed in RDKit's
/// order so the value is bit-identical to `Crippen.MolLogP`.
pub fn logp_crippen(mol: &Molecule) -> f64 {
    crippen_totals(mol).0
}

// ---------------------------------------------------------------------------
// 9. Lipinski Rule of Five
// ---------------------------------------------------------------------------

/// Apply Lipinski's Rule of Five.
///
/// Returns `true` when all four criteria are satisfied:
/// - Molecular weight ≤ 500 Da
/// - H-bond donors ≤ 5
/// - H-bond acceptors ≤ 10
/// - Crippen LogP ≤ 5.0
pub fn lipinski_passes(mol: &Molecule) -> bool {
    molecular_weight(mol) <= 500.0
        && hbd_count(mol) <= 5
        && hba_count(mol) <= 10
        && logp_crippen(mol) <= 5.0
}

// ---------------------------------------------------------------------------
// Fsp3 — fraction of sp3 carbons
// ---------------------------------------------------------------------------

/// Fraction of sp3 carbons: sp3_C / total_C.
///
/// sp3 carbon is defined as a neutral, non-aromatic carbon with a complete
/// valence of four and no double or triple bond to any neighbour. The
/// valence check is important for bracket atoms such as [13C]: RDKit does
/// not count zero-valence carbon radicals as CSP3.
/// Returns 0.0 if the molecule contains no carbon atoms.
pub fn fsp3(mol: &Molecule) -> f64 {
    chematic_smiles::with_rdkit_model_molecule(mol, fsp3_impl)
}

fn fsp3_impl(mol: &Molecule) -> f64 {
    // RDKit `calcFractionCSP3`: carbons whose total degree (neighbours plus
    // hydrogens) is four, over all carbons. A four-bonded carbanion (the
    // ring carbon of a ferrocenyl `[C-]` bonded to Fe) counts.
    let is_c = |a: &chematic_core::Atom| !a.wildcard && a.element.atomic_number() == 6;
    let c_total = mol.atoms().filter(|(_, a)| is_c(a)).count();
    if c_total == 0 {
        return 0.0;
    }
    let sp3 = mol
        .atoms()
        .filter(|(idx, a)| {
            is_c(a) && mol.degree(*idx) + usize::from(mol.implicit_hydrogen_count(*idx)) == 4
        })
        .count();
    sp3 as f64 / c_total as f64
}

// ---------------------------------------------------------------------------
// Aromatic ring count
// ---------------------------------------------------------------------------

/// Number of aromatic rings in the molecule (from SSSR).
///
/// A ring is considered aromatic when every atom in it carries the
/// `aromatic` flag.
pub fn aromatic_ring_count(mol: &Molecule) -> usize {
    use chematic_perception::count_aromatic_rings;
    count_aromatic_rings(mol)
}

/// Count aromatic rings after applying the opt-in RDKit aromaticity model.
/// Native aromatic flags and native ring counts are left untouched.
pub fn rdkit_aromatic_ring_count(mol: &Molecule) -> usize {
    if chematic_perception::rdkit_parity_view_is_identity(mol) {
        return chematic_perception::aromatic_ring_list_preperceived(mol).len();
    }
    match chematic_perception::apply_aromaticity_rdkit_parity_shared(mol).as_ref() {
        // The descriptor view reads aromatic-flagged ring triple bonds as
        // aromatic, as RDKit's ring counting does.
        Ok(_) => {
            chematic_perception::aromatic_ring_list_preperceived(&descriptor_aromaticity(mol)).len()
        }
        Err(_) => chematic_perception::aromatic_ring_list(mol).len(),
    }
}

// ---------------------------------------------------------------------------
// 11. Formal charge sum
// ---------------------------------------------------------------------------

/// Sum of all formal charges in the molecule.
pub fn formal_charge_sum(mol: &Molecule) -> i32 {
    mol.atoms().map(|(_, a)| a.charge as i32).sum()
}

// ---------------------------------------------------------------------------
// 12. Molar Refractivity (Wildman-Crippen additive model)

/// Return the MR contribution per implicit-H attached to `parent_idx`.
fn h_mr_for_parent(
    mol: &Molecule,
    parent_idx: AtomIdx,
    parent_an: u8,
    parent_aromatic: bool,
    fallback: f64,
) -> f64 {
    match parent_an {
        6 => 1.057,  // H1: H on C
        7 => 0.9627, // H3: H on N
        8 => {
            if parent_aromatic {
                fallback // aromatic O — HS fallback
            } else if oh_hydrogen_is_h2(mol, parent_idx) {
                1.395
            } else if mol.neighbors(parent_idx).any(|(nb, _)| {
                // [#1]O[#7]: N-O-H (hydroxamate, oxime) → same MR as H on N
                mol.atom(nb).element.atomic_number() == 7
            }) {
                0.9627
            } else if mol.neighbors(parent_idx).any(|(nb, _)| {
                // [#1]OC=[#6,#7,O,S]: enol, carboxylic acid, or O[O,S]
                let an = mol.atom(nb).element.atomic_number();
                (an == 6
                    && mol.neighbors(nb).any(|(nb2, bidx2)| {
                        nb2 != parent_idx
                            && mol.bond(bidx2).order == BondOrder::Double
                            && matches!(mol.atom(nb2).element.atomic_number(), 6 | 7 | 8 | 16)
                    }))
                    || an == 8
                    || an == 16
            }) {
                1.805
            } else {
                1.395 // H2: aliphatic or phenolic OH
            }
        }
        _ => 1.395, // [#1][!C;!N;!O] → 1.395 (SH, PH, etc.)
    }
}

/// Per-atom Molar Refractivity contributions (Wildman & Crippen 1999, RDKit-compatible).
///
/// Uses the same 117-entry SMARTS atom-type table as `logp_crippen_per_atom` but
/// reads the MR column. Results match RDKit `Crippen.MolMR()` for common molecules.
/// H contributions are folded into the attached heavy atom. Index matches mol.atoms().
pub fn mr_per_atom(mol: &Molecule) -> Vec<f64> {
    let queries = get_crippen_queries();
    let mol_arom = descriptor_aromaticity(mol);
    let anchor_types = crippen_anchor_types(&mol_arom, queries);

    let h_fallback = CRIPPEN_SMARTS
        .iter()
        .find(|(sma, _, _)| *sma == "[#1]")
        .map(|e| e.2)
        .unwrap_or(1.112);

    mol.atoms()
        .map(|(idx, atom)| {
            if atom.element.atomic_number() == 1 {
                return 0.0;
            }
            let h_count_impl = implicit_hcount(mol, idx);
            let h_count_expl = mol
                .neighbors(idx)
                .filter(|(nb, _)| mol.atom(*nb).element.atomic_number() == 1)
                .count() as u8;
            let h_count = h_count_impl + h_count_expl;
            let heavy = anchor_types[idx.0 as usize]
                .map(|t| queries[t].2)
                .unwrap_or(0.0);
            // RDKit adds each hydrogen's contribution to its heavy atom one
            // at a time (getCrippenAtomContribs folds the added H atoms
            // back); a single multiplication rounds differently.
            let mut contrib = heavy;
            if h_count > 0 {
                let h = h_mr_for_parent(
                    &mol_arom,
                    idx,
                    atom.element.atomic_number(),
                    // Hydrogen typing sees the RDKit-perceived environment, like
                    // the heavy-atom SMARTS typing (Kekule phenol OH is H2, not H4).
                    mol_arom.atom(idx).aromatic,
                    h_fallback,
                );
                for _ in 0..h_count {
                    contrib += h;
                }
            }
            contrib
        })
        .collect()
}

/// Compute Molar Refractivity using the Wildman-Crippen additive model.
///
/// Uses the same atom-type framework as `logp_crippen` but with MR contributions
/// from Wildman & Crippen 1999 (J. Chem. Inf. Comput. Sci. 39, 868-873).
pub fn molar_refractivity(mol: &Molecule) -> f64 {
    crippen_totals(mol).1
}

/// Compute both LogP and MR from a single Crippen anchor-set pass.
///
/// Equivalent to calling `logp_crippen(mol)` and `molar_refractivity(mol)` in
/// sequence but shares the 117-pattern SMARTS matching (including the single SSSR
/// computation), making it roughly 2× faster when both values are needed.
pub fn logp_and_mr(mol: &Molecule) -> (f64, f64) {
    crippen_totals(mol)
}

/// Crippen LogP and MR totals, summed in the order RDKit's
/// `calcCrippenDescriptors` sums them, so the totals are bit-identical:
/// RDKit adds hydrogens (`AddHs` appends them after every graph atom, atom
/// by atom) and sums every atom's contribution in index order. So the graph
/// atoms come first (an explicit hydrogen typed by its neighbour), then each
/// heavy atom's implicit hydrogens in heavy-atom order.
fn crippen_totals(mol: &Molecule) -> (f64, f64) {
    // Molecules chematic's perception may model differently from RDKit
    // (radicals, dummies, hydrogen graph atoms, unusual valences, ...):
    // type RDKit's own `AddHs(MolFromSmiles(s))` molecule with the same
    // table, as `Crippen.MolLogP` does.
    if crippen_needs_rdkit_model(mol)
        && let Some(t) = crippen_totals_rdkit_model(mol)
    {
        return t;
    }
    let queries = get_crippen_queries();
    let mol_arom = descriptor_aromaticity(mol);
    let anchor_types = crippen_anchor_types(&mol_arom, queries);
    let (h_logp_fallback, h_mr_fallback) = CRIPPEN_SMARTS
        .iter()
        .find(|(sma, _, _)| *sma == "[#1]")
        .map(|e| (e.1, e.2))
        .unwrap_or((0.1125, 1.112));
    // Contribution of one hydrogen on `parent`.
    let h_on = |parent: AtomIdx| {
        let an = mol.atom(parent).element.atomic_number();
        // Hydrogen typing sees the RDKit-perceived environment, like the
        // heavy-atom SMARTS typing (Kekule phenol OH is H2, not H4).
        let aromatic = mol_arom.atom(parent).aromatic;
        (
            h_logp_for_parent(&mol_arom, parent, an, aromatic, h_logp_fallback),
            h_mr_for_parent(&mol_arom, parent, an, aromatic, h_mr_fallback),
        )
    };

    let (mut logp, mut mr) = (0.0f64, 0.0f64);
    for (idx, atom) in mol.atoms() {
        if atom.element.atomic_number() == 1 {
            // A graph hydrogen keeps its own slot (0.0 when unattached, as
            // the per-atom API reports it).
            if let Some((parent, _)) = mol
                .neighbors(idx)
                .find(|(nb, _)| mol.atom(*nb).element.atomic_number() != 1)
            {
                let (l, m) = h_on(parent);
                logp += l;
                mr += m;
            }
            continue;
        }
        logp += anchor_types[idx.0 as usize]
            .map(|t| queries[t].1)
            .unwrap_or(0.0);
        mr += anchor_types[idx.0 as usize]
            .map(|t| queries[t].2)
            .unwrap_or(0.0);
    }
    for (idx, atom) in mol.atoms() {
        if atom.element.atomic_number() == 1 {
            continue;
        }
        let n = implicit_hcount(mol, idx);
        if n > 0 {
            let (l, m) = h_on(idx);
            for _ in 0..n {
                logp += l;
                mr += m;
            }
        }
    }
    (logp, mr)
}
