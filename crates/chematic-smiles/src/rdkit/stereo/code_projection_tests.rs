//! Frozen baseline comparison; only `_ChiralityPossible` may differ in the projection.
use super::super::{Molecule, RdkitLegacyStereo, RdkitSmilesError, mol, parse, sanitize};
use super::*;

fn compare_helpers(m: &Mol) {
    let n = m.atoms.len();
    if n == 0 {
        return;
    }
    for ranks in [
        vec![0; n],
        (0..n as u32).collect(),
        (0..n).map(|i| (i % 3) as u32 % n as u32).collect(),
    ] {
        for a in 0..n {
            let mut old = vec![(123, usize::MAX)];
            let mut new = old.clone();
            assert_eq!(
                frozen_potential(m, a, &ranks, &mut old),
                is_atom_potential_chiral_center(m, a, &ranks, &mut new)
            );
            assert_eq!(
                old, new,
                "neighbor append order/duplicate break at atom {a}"
            );
        }
    }
}

fn compare_projection(m: &Mol) {
    let mut old = m.clone();
    let mut new = m.clone();
    frozen_legacy(&mut old, true, true, false);
    legacy_stereo_perception_unflagged(&mut new);
    for a in &mut old.atoms {
        a.chirality_possible = false;
    }
    for a in &mut new.atoms {
        a.chirality_possible = false;
    }
    assert_eq!(format!("{old:?}"), format!("{new:?}"));
}

#[test]
fn full_corpus_projection_and_neighbor_order_match_frozen_baseline() {
    let corpora = [
        include_str!("../../../../../validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi"),
        include_str!("../../../../../validation/cosmolkit_comparison/corpora/stress-v1.smi"),
    ];
    let mut rows = 0;
    let mut raw_views = 0;
    let mut sanitized_views = 0;
    for line in corpora
        .iter()
        .flat_map(|s| s.lines())
        .filter(|s| !s.trim().is_empty())
    {
        let text = line.split_whitespace().next().unwrap();
        let core = crate::parse(text).unwrap();
        rows += 1;
        assert_eq!(
            format!("{:?}", frozen_api(&core)),
            format!("{:?}", super::super::rdkit_legacy_stereo(&core)),
            "API input {text}"
        );
        let Ok(mut raw) = parse::from_chematic(&core) else {
            continue;
        };
        let mut sanitized = raw.clone();
        if raw.find_rings().is_ok() {
            compare_helpers(&raw);
            compare_projection(&raw);
            raw_views += 1;
        }
        if sanitize::remove_hs_and_sanitize(&mut sanitized).is_ok() {
            compare_helpers(&sanitized);
            compare_projection(&sanitized);
            sanitized_views += 1;
        }
    }
    assert_eq!(rows, 14_072);
    assert!(raw_views >= 10_000 && sanitized_views >= 10_000);
    eprintln!("frozen API rows={rows}; raw views={raw_views}; sanitized views={sanitized_views}");
}

#[test]
fn dative_high_degree_and_full_flag_modes_keep_frozen_results() {
    let mut cases = 0;
    for anum in [6, 7, 15, 33, 16, 34] {
        for covalent in [2, 3, 4, 5] {
            for dative in [0, 1, 8, 64] {
                for outgoing in [false, true] {
                    for hs in [0, 1] {
                        let mut m = Mol::default();
                        let mut atom = mol::Atom::new(anum);
                        atom.chiral = ChiralTag::Ccw;
                        atom.num_explicit_hs = hs;
                        atom.no_implicit = true;
                        atom.explicit_valence = covalent as i32 + hs as i32;
                        atom.hybrid = Hybridization::Sp3;
                        m.add_atom(atom);
                        for i in 0..covalent + dative {
                            let mut atom = mol::Atom::new([6, 7, 8, 9][i % 4]);
                            atom.no_implicit = true;
                            let a = m.add_atom(atom);
                            let (begin, end) = if outgoing { (0, a) } else { (a, 0) };
                            m.add_bond(mol::Bond::new(
                                begin,
                                end,
                                if i < covalent {
                                    BondType::Single
                                } else {
                                    BondType::Dative
                                },
                            ));
                        }
                        m.find_rings().unwrap();
                        compare_helpers(&m);
                        compare_projection(&m);
                        for clean in [false, true] {
                            for flag in [false, true] {
                                let mut old = m.clone();
                                let mut new = m.clone();
                                assert_eq!(
                                    frozen_legacy(&mut old, clean, flag, false),
                                    legacy_stereo_perception(&mut new, clean, flag)
                                );
                                assert_eq!(format!("{old:?}"), format!("{new:?}"));
                            }
                        }
                        cases += 1;
                    }
                }
            }
        }
    }
    assert_eq!(cases, 384);
}

#[test]
fn empty_projection_matches_frozen_baseline() {
    let mut m = Mol::default();
    m.find_rings().unwrap();
    compare_projection(&m);
}

// Frozen functions below: only names and formatting differ.

fn frozen_potential(
    mol: &Mol,
    a: usize,
    ranks: &[u32],
    nbrs: &mut Vec<(u32, usize)>,
) -> (bool, bool) {
    let atom = &mol.atoms[a];
    let mut legal = true;
    let mut dupes = false;
    let nz = atom_nonzero_degree(mol, a);
    let tnz = nz + mol.total_num_hs(a) as usize;
    if tnz > 4 {
        legal = false;
    } else {
        if tnz < 3 {
            legal = false;
        } else if nz < 3 && atom.anum != 15 && atom.anum != 33 {
            legal = false;
        } else if nz == 3 {
            if mol.total_num_hs(a) == 1 {
                if has_protium_neighbor(mol, a) {
                    legal = false;
                }
            } else {
                legal = false;
                if atom.anum == 7 {
                    if atom.hybrid == Hybridization::Sp3
                        && !atom_has_conjugated_bond(mol, a)
                        && (mol.ring_info().is_atom_in_ring_of_size(a, 3)
                            || is_atom_bridgehead(mol, a))
                    {
                        legal = true;
                    }
                } else if atom.anum == 15 || atom.anum == 33 {
                    legal = true;
                } else if (atom.anum == 16 || atom.anum == 34)
                    && (atom.explicit_valence == 4
                        || (atom.explicit_valence == 3 && atom.charge == 1))
                {
                    legal = true;
                }
            }
        }
        if legal && !ranks.is_empty() {
            let mut seen = vec![false; mol.atoms.len()];
            for &b in &mol.atom_bonds[a] {
                let other = mol.bonds[b].other(a);
                nbrs.push((ranks[other], b));
                if !bond_affects_atom_chirality(mol, b, a) {
                    continue;
                }
                let r = ranks[other] as usize;
                if seen[r] {
                    dupes = true;
                    break;
                }
                seen[r] = true;
            }
        }
    }
    (legal, dupes)
}

fn frozen_assign(mol: &mut Mol, ranks: &mut Vec<u32>, flag_possible: bool) -> (bool, bool) {
    let mut atom_changed = false;
    let mut unassigned = 0;
    for a in 0..mol.atoms.len() {
        let tag = mol.atoms[a].chiral;
        if !(flag_possible || tag != ChiralTag::Unspecified) {
            continue;
        }
        if mol.atoms[a].cip_code.is_some() {
            continue;
        }
        if ranks.is_empty() {
            assign_atom_cip_ranks(mol, ranks);
        }
        let mut nbrs = Vec::new();
        let (legal, dupes) = frozen_potential(mol, a, ranks, &mut nbrs);
        if legal {
            unassigned += 1;
        }
        if legal && !dupes && flag_possible {
            mol.atoms[a].chirality_possible = true;
        }
        if legal && !dupes && tag != ChiralTag::Unspecified {
            atom_changed = true;
            unassigned -= 1;
            // std::sort with Rankers::pairLess (ranks are distinct here,
            // except for bonds that do not affect chirality).
            nbrs.sort_by_key(|p| p.0);
            let indices: Vec<usize> = nbrs.iter().map(|p| p.1).collect();
            let mut odd = mol.perturbation_is_odd(a, &indices);
            if indices.len() == 3 && mol.total_num_hs(a) == 1 {
                odd = !odd;
            }
            let mut t = tag;
            if odd {
                t = if t == ChiralTag::Ccw {
                    ChiralTag::Cw
                } else {
                    ChiralTag::Ccw
                };
            }
            mol.atoms[a].cip_code = Some(if t == ChiralTag::Ccw { b'S' } else { b'R' });
        }
    }
    (unassigned > 0, atom_changed)
}

fn frozen_legacy(
    mol: &mut Mol,
    clean_it: bool,
    flag_possible: bool,
    skip_flag_only_pass: bool,
) -> Vec<u32> {
    let mut has_stereo_atoms = false;
    let mut has_potential_stereo_atoms = false;
    for a in 0..mol.atoms.len() {
        if clean_it {
            mol.atoms[a].cip_code = None;
            mol.atoms[a].chirality_possible = false;
        }
        if !has_stereo_atoms && mol.atoms[a].chiral != ChiralTag::Unspecified {
            has_stereo_atoms = true;
        } else if !has_potential_stereo_atoms {
            let mut nbrs = Vec::new();
            has_potential_stereo_atoms = frozen_potential(mol, a, &[], &mut nbrs).0;
        }
    }
    let mut has_stereo_bonds = false;
    let mut has_potential_stereo_bonds = false;
    for b in 0..mol.bonds.len() {
        if clean_it {
            let bt = mol.bonds[b].bt;
            if (bt == BondType::Double || bt == BondType::Aromatic)
                && !should_detect_double_bond_stereo(mol, b)
            {
                if mol.bonds[b].stereo != BondStereo::None {
                    mol.bonds[b].stereo = BondStereo::None;
                    mol.bonds[b].stereo_atoms.clear();
                }
                continue;
            } else if bt == BondType::Double && mol.bonds[b].stereo != BondStereo::Any {
                mol.bonds[b].stereo = BondStereo::None;
                mol.bonds[b].stereo_atoms.clear();
            }
        }
        if !has_stereo_bonds && mol.bonds[b].bt == BondType::Double {
            let mut is_specified = false;
            if mol.bonds[b].requested.is_some() {
                has_stereo_bonds = true;
                is_specified = true;
            }
            for end in [mol.bonds[b].begin, mol.bonds[b].end] {
                if mol.atom_bonds[end]
                    .iter()
                    .any(|&nb| mol.bonds[nb].dir.is_set())
                {
                    has_stereo_bonds = true;
                    is_specified = true;
                    break;
                }
            }
            if !has_potential_stereo_bonds
                && !is_specified
                && should_detect_double_bond_stereo(mol, b)
            {
                has_potential_stereo_bonds = true;
            }
        }
        if !clean_it && has_stereo_bonds && has_potential_stereo_bonds {
            break;
        }
    }
    let mut ranks: Vec<u32> = Vec::new();
    let mut keep_going = has_stereo_atoms | has_stereo_bonds;
    if !keep_going && !skip_flag_only_pass {
        keep_going = flag_possible && (has_potential_stereo_atoms || has_potential_stereo_bonds);
    }
    while keep_going {
        let changed_atoms;
        let changed_bonds;
        if has_stereo_atoms || has_potential_stereo_atoms {
            let (h, c) = frozen_assign(mol, &mut ranks, flag_possible);
            has_stereo_atoms = h;
            changed_atoms = c;
        } else {
            changed_atoms = false;
        }
        if has_stereo_bonds || has_potential_stereo_bonds {
            let (h, c) = assign_bond_stereo_codes(mol, &mut ranks);
            has_stereo_bonds = h;
            changed_bonds = c;
        } else {
            changed_bonds = false;
        }
        keep_going = (has_stereo_atoms || has_stereo_bonds) && (changed_atoms || changed_bonds);
        if keep_going {
            rerank_atoms(mol, &mut ranks);
        }
    }
    if !clean_it {
        return ranks;
    }
    for atom in &mut mol.atoms {
        atom.ring_stereochem_cand = None;
        atom.ring_stereo_atoms = None;
    }
    let possible =
        if ranks.is_empty() && !mol.atoms.iter().any(|a| a.chiral != ChiralTag::Unspecified) {
            vec![false; mol.atoms.len()]
        } else {
            if ranks.is_empty() {
                assign_atom_cip_ranks(mol, &mut ranks);
            }
            find_chiral_atom_special_cases(mol, &ranks)
        };
    for a in 0..mol.atoms.len() {
        let atom = &mol.atoms[a];
        if atom.chiral != ChiralTag::Unspecified
            && atom.chiral.nontet().is_none()
            && atom.cip_code.is_none()
            && (!possible[a] || atom.ring_stereo_atoms.is_none())
        {
            mol.atoms[a].chiral = ChiralTag::Unspecified;
            let atom = &mol.atoms[a];
            if atom.num_explicit_hs == 1 && atom.charge == 0 && !atom.aromatic {
                mol.atoms[a].num_explicit_hs = 0;
                mol.atoms[a].no_implicit = false;
                // calcExplicitValence(false) / calcImplicitValence(false)
                if let Ok(v) = mol.calc_explicit_valence_value(a, false) {
                    mol.atoms[a].explicit_valence = v;
                }
                if let Ok(v) = mol.calc_implicit_valence_value(a, false) {
                    mol.atoms[a].implicit_valence = v;
                }
            }
        }
    }
    for b in 0..mol.bonds.len() {
        let (beg, end) = (mol.bonds[b].begin, mol.bonds[b].end);
        if mol.bonds[b].bt == BondType::Double
            && mol.bonds[b].stereo == BondStereo::Any
            && (mol.degree(beg) == 1 || mol.degree(end) == 1)
        {
            mol.bonds[b].stereo = BondStereo::None;
        }
        if mol.bonds[b].bt == BondType::Double
            && matches!(mol.bonds[b].stereo, BondStereo::Any | BondStereo::None)
        {
            for batom in [beg, end] {
                for nb_i in mol.atom_bonds[batom].clone() {
                    if nb_i == b {
                        continue;
                    }
                    let nbond = &mol.bonds[nb_i];
                    if nbond.dir.is_set()
                        && matches!(nbond.bt, BondType::Single | BondType::Aromatic)
                    {
                        let far = nbond.other(batom);
                        let ok_to_clear = !mol.atom_bonds[far].iter().any(|&nb_j| {
                            mol.bonds[nb_j].bt == BondType::Double
                                && !matches!(
                                    mol.bonds[nb_j].stereo,
                                    BondStereo::Any | BondStereo::None
                                )
                        });
                        if ok_to_clear {
                            mol.bonds[nb_i].dir = BondDir::None;
                        }
                    }
                }
            }
        }
    }
    ranks
}

fn frozen_api(mol: &Molecule) -> Result<RdkitLegacyStereo, RdkitSmilesError> {
    let mut m = parse::from_chematic(mol)?;
    let n_atoms = m.atoms.len();
    sanitize::remove_hs_and_sanitize(&mut m)?;
    if m.atoms.len() != n_atoms || m.bonds.len() != mol.bond_count() {
        return Err(RdkitSmilesError::Unsupported(
            "hydrogen removal renumbers atoms".into(),
        ));
    }
    frozen_legacy(&mut m, true, true, false);
    let order = mol.rdkit_bond_order();
    let mut bond_stereo = vec![0u8; mol.bond_count()];
    for (k, b) in order.iter().enumerate() {
        bond_stereo[b.0 as usize] = m.bonds[k].stereo as u8;
    }
    Ok(RdkitLegacyStereo {
        atom_tagged: m
            .atoms
            .iter()
            .map(|a| a.chiral != mol::ChiralTag::Unspecified)
            .collect(),
        atom_cip: m.atoms.iter().map(|a| a.cip_code).collect(),
        bond_stereo,
    })
}
