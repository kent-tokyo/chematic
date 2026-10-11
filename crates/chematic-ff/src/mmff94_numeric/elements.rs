// ── C type assignment ────────────────────────────────────────────────────────

fn assign_c_type(
    mol: &Molecule,
    rings: &[Vec<AtomIdx>],
    idx: AtomIdx,
) -> Result<u8, NumericTypeError> {
    let atom = mol.atom(idx);

    // Aromatic: detect 5-ring vs 6-ring (RDKit AtomTyper.cpp, aromatic block)
    if atom.aromatic
        && let Some(t) = aromatic_c_type(mol, rings, idx)
    {
        return Ok(t);
    }
    // Aromatic-flagged but not in a fully-aromatic-bonded 5- or 6-ring (e.g.
    // a Kekule-only representation slipping through, or an aromatic ring
    // size MMFF94 doesn't special-case) -- fall through to the aliphatic
    // rules below rather than silently mis-assign a wrong-element type,
    // matching RDKit's own fallthrough when its aromatic switch doesn't set
    // atomType.

    // Aliphatic carbon, dispatched on total degree exactly as RDKit's
    // `setMMFFHeavyAtomType` (`AtomTyper.cpp`, Release_2026_03_6, `case 6`
    // aliphatic block). The previous port keyed on "has a double bond" and
    // fell back to alkyl carbon (1), so a three-connected carbon with no
    // double bond -- a ring carbon of an MMFF-aromatic ring RDKit's typer
    // does not special-case (tropone's seven-ring) -- was typed 1 instead of
    // 2, and CR4E (30), CO2M/CS2M (41) and CNN+/CGD+ (57) were never
    // assigned.
    match total_degree(mol, idx) {
        4 => {
            if atom_in_ring_of_size(rings, idx, 3) {
                return Ok(22); // CR3R
            }
            if atom_in_ring_of_size(rings, idx, 4) {
                return Ok(20); // CR4R
            }
            Ok(1) // CR alkyl carbon
        }
        3 => {
            let (mut n_n2, mut n_n3, mut n_o, mut n_s) = (0, 0, 0, 0);
            let mut double_bonded_element: Option<Element> = None;
            for nb in bonds_of(mol, idx) {
                let nb_atom = mol.atom(nb.neighbor);
                if nb.order == BondOrder::Double {
                    double_bonded_element = Some(nb_atom.element);
                }
                let nb_degree = total_degree(mol, nb.neighbor);
                if nb_degree == 1 {
                    match nb_atom.element {
                        Element::O => n_o += 1,
                        Element::S => n_s += 1,
                        _ => {}
                    }
                } else if nb_atom.element == Element::N {
                    if nb_degree == 3 {
                        n_n3 += 1;
                    } else if nb_degree == 2 && nb.order == BondOrder::Double {
                        n_n2 += 1;
                    }
                }
            }
            if n_n3 >= 2 && n_n2 == 0 && double_bonded_element == Some(Element::N) {
                return Ok(57); // CNN+ / CGD+
            }
            if n_o == 2 || n_s == 2 {
                return Ok(41); // CO2M / CS2M
            }
            if atom_in_ring_of_size(rings, idx, 4) && double_bonded_element == Some(Element::C) {
                return Ok(30); // CR4E
            }
            if matches!(
                double_bonded_element,
                Some(Element::N | Element::O | Element::P | Element::S)
            ) {
                return Ok(3); // C=N / C=O / C=P / C=S family
            }
            Ok(2) // C=C vinylic / generic sp2 carbon
        }
        2 => Ok(4), // CSP: acetylenic or allenic carbon
        _ => Ok(1),
    }
}

/// Faithful port of RDKit's aromatic-carbon cases (`AtomTyper.cpp`
/// `setMMFFHeavyAtomType`, `case 6:` under both the 5-ring and 6-ring
/// aromatic switches). Returns `None` if the atom is flagged aromatic but
/// isn't actually in a fully-aromatic-bonded 5- or 6-membered ring.
///
/// CIM+ (type 80): with no beta heteroatom, a carbon bonded to two
/// three-connected N (one in an aromatic 5-ring, none in an aromatic 6-ring)
/// of which one is a positive non-N-oxide N, or to three such N with two in
/// aromatic 5-rings.
fn aromatic_c_type(mol: &Molecule, rings: &[Vec<AtomIdx>], idx: AtomIdx) -> Option<u8> {
    if atom_in_aromatic_ring_of_size(mol, rings, idx, 5)
        && let Some(atom_type) = aromatic_five_ring_c_type(mol, rings, idx)
    {
        return Some(atom_type);
    }

    if atom_in_aromatic_ring_of_size(mol, rings, idx, 6) {
        return Some(37); // CB: benzene/pyridine-ring carbon
    }

    None
}

fn aromatic_five_ring_c_type(
    mol: &Molecule,
    rings: &[Vec<AtomIdx>],
    idx: AtomIdx,
) -> Option<u8> {
    let hetero = find_alpha_beta_heteroatoms(mol, rings, idx);
    if hetero.beta.is_empty() && is_cim_plus_carbon(mol, rings, idx) {
        return Some(80);
    }
    if hetero.alpha.len() == hetero.beta.len()
        && is_general_c5_carbon(mol, rings, idx, &hetero)
    {
        return Some(78);
    }
    if !hetero.alpha.is_empty() && (hetero.beta.is_empty() || hetero.is_alpha_os) {
        return Some(63);
    }
    if !hetero.beta.is_empty() && (hetero.alpha.is_empty() || hetero.is_beta_os) {
        return Some(64);
    }
    None
}

fn is_cim_plus_carbon(mol: &Molecule, rings: &[Vec<AtomIdx>], idx: AtomIdx) -> bool {
    let mut counts = [0usize; 4]; // N, charged N, aromatic-5 N, aromatic-6 N
    for neighbor in bonds_of(mol, idx) {
        let atom = neighbor.neighbor;
        if mol.atom(atom).element != Element::N || total_degree(mol, atom) != 3 {
            continue;
        }
        counts[0] += 1;
        counts[1] += usize::from(mol.atom(atom).charge > 0 && !is_atom_n_oxide(mol, atom));
        counts[2] += usize::from(atom_in_aromatic_ring_of_size(mol, rings, atom, 5));
        counts[3] += usize::from(atom_in_aromatic_ring_of_size(mol, rings, atom, 6));
    }
    ((counts[0] == 2 && counts[2] > 0) || (counts[0] == 3 && counts[2] == 2))
        && counts[1] > 0
        && counts[3] == 0
}

fn is_general_c5_carbon(
    mol: &Molecule,
    rings: &[Vec<AtomIdx>],
    idx: AtomIdx,
    hetero: &AlphaBetaHeteroatoms,
) -> bool {
    let neighbors = bonds_of(mol, idx);
    let benzene_surrounded = neighbors.iter().all(|neighbor| {
        mol.atom(neighbor.neighbor).element == Element::C
            && atom_in_ring_of_size(rings, neighbor.neighbor, 6)
    });
    let aromatic_surrounded = neighbors.iter().all(|neighbor| {
        !rings.iter().any(|ring| {
            ring.len() == 5 && ring.contains(&idx) && ring.contains(&neighbor.neighbor)
        }) || mol.atom(neighbor.neighbor).aromatic
    });
    (hetero.alpha.is_empty()
        && hetero.beta.is_empty()
        && !benzene_surrounded
        && aromatic_surrounded)
        || (!hetero.alpha.is_empty()
            && !hetero.beta.is_empty()
            && (!hetero.alpha_or_beta_in_same_ring
                || (!hetero.is_alpha_os && !hetero.is_beta_os)))
}

// ── N type assignment ────────────────────────────────────────────────────────

fn assign_n_type(
    mol: &Molecule,
    rings: &[Vec<AtomIdx>],
    idx: AtomIdx,
) -> Result<u8, NumericTypeError> {
    let atom = mol.atom(idx);

    // Aromatic nitrogen
    if atom.aromatic
        && let Some(t) = aromatic_n_type(mol, rings, idx)
    {
        return Ok(t);
    }

    Ok(general_n_type(mol, rings, idx))
}

/// RDKit's `getValence(EXPLICIT) + getNumImplicitHs()` for a nitrogen of the
/// MMFF view. `setMMFFAromaticity` sets the bonds of an accepted ring to
/// `AROMATIC` and recomputes the explicit valence of its non-carbon atoms
/// from them (`calculateExplicitValence`: aromatic bonds count 1.5, an
/// aromatic atom over its default valence drops to the nearest allowed
/// valence within 1.5, then `round(accum + 0.1)`); every other atom keeps
/// the valence of RDKit's Kekule structure.
fn rdkit_n_total_bond_order(mol: &Molecule, idx: AtomIdx) -> u32 {
    let atom = mol.atom(idx);
    let nbrs = bonds_of(mol, idx);
    let implicit = u32::from(implicit_hcount(mol, idx));
    if !atom.aromatic || !nbrs.iter().any(|b| b.order == BondOrder::Aromatic) {
        return nbrs.iter().map(|b| b.order.order_int() as u32).sum::<u32>() + implicit;
    }
    let mut accum: f64 = nbrs
        .iter()
        .map(|b| match b.order {
            BondOrder::Aromatic => 1.5,
            o => f64::from(o.order_int()),
        })
        .sum::<f64>()
        + f64::from(implicit);
    // Effective element: N+ is valence-isoelectronic with C, N- with O.
    let (dv, valens): (f64, &[f64]) = match atom.charge {
        1 => (4.0, &[4.0]),
        -1 => (2.0, &[2.0]),
        0 => (3.0, &[3.0]),
        _ => (3.0, &[3.0]),
    };
    if accum > dv {
        let mut pval = dv;
        for &v in valens {
            if v > accum {
                break;
            }
            pval = v;
        }
        if accum - pval <= 1.5 {
            accum = pval;
        }
    }
    (accum + 0.1).round() as u32
}

/// Literal port of RDKit's non-aromatic nitrogen typing (`AtomTyper.cpp`
/// `setMMFFHeavyAtomType`, `case 7:` of the general switch, RDKit
/// 2026.03.6), including its loop-carried state: `elementDoubleBondedToC`,
/// `elementTripleBondedToC` and `isNbrBenzeneC` keep their values from one
/// carbon neighbour to the next, so a carbon met later can be judged with an
/// earlier carbon's double bond (a neutral N between an imine carbon and a
/// tri-aminated carbon is NGD+, 56, in RDKit). Neighbours are visited in
/// bond order, which is RDKit's order for molecules read from MOL files and
/// for SMILES without ring-closure bonds on the visited atoms.
fn general_n_type(mol: &Molecule, rings: &[Vec<AtomIdx>], idx: AtomIdx) -> u8 {
    let neighbors = bonds_of(mol, idx);
    let degree = total_degree(mol, idx);
    let total_bond_order = rdkit_n_total_bond_order(mol, idx);
    let (terminal_oxygen, mut nso_family) =
        nitrogen_terminal_environment(mol, &neighbors, total_bond_order);

    if degree == 4 {
        return if is_atom_n_oxide(mol, idx) { 68 } else { 34 };
    }
    if degree == 3
        && let Some(atom_type) = degree_three_n_type(
            mol,
            rings,
            idx,
            &neighbors,
            total_bond_order,
            terminal_oxygen,
            &mut nso_family,
        )
    {
        return atom_type;
    }
    if degree == 2
        && let Some(atom_type) = degree_two_n_type(
            mol,
            idx,
            &neighbors,
            total_bond_order,
            terminal_oxygen,
            nso_family,
        )
    {
        return atom_type;
    }
    if nso_family {
        return 43;
    }
    if degree == 1 && let Some(atom_type) = degree_one_n_type(mol, idx, &neighbors) {
        return atom_type;
    }
    8
}

fn nitrogen_terminal_environment(
    mol: &Molecule,
    neighbors: &[BondInfo],
    total_bond_order: u32,
) -> (usize, bool) {
    let mut terminal_oxygen = 0;
    let mut nso_family = false;
    for neighbor in neighbors {
        let element = mol.atom(neighbor.neighbor).element;
        terminal_oxygen += usize::from(
            element == Element::O && total_degree(mol, neighbor.neighbor) == 1,
        );
        if total_bond_order >= 3 && matches!(element, Element::P | Element::S) {
            nso_family |= bonds_of(mol, neighbor.neighbor)
                .iter()
                .filter(|bond| {
                    mol.atom(bond.neighbor).element == Element::O
                        && total_degree(mol, bond.neighbor) == 1
                })
                .count()
                >= 2;
        }
    }
    (terminal_oxygen, nso_family)
}

fn degree_three_n_type(
    mol: &Molecule,
    rings: &[Vec<AtomIdx>],
    idx: AtomIdx,
    neighbors: &[BondInfo],
    total_bond_order: u32,
    terminal_oxygen: usize,
    nso_family: &mut bool,
) -> Option<u8> {
    if total_bond_order >= 4
        && let Some(atom_type) =
            charged_degree_three_n_type(mol, idx, neighbors, terminal_oxygen)
    {
        return Some(atom_type);
    }
    (total_bond_order >= 3)
        .then(|| conjugated_degree_three_n_type(mol, rings, neighbors, nso_family))
        .flatten()
}

fn charged_degree_three_n_type(
    mol: &Molecule,
    idx: AtomIdx,
    neighbors: &[BondInfo],
    terminal_oxygen: usize,
) -> Option<u8> {
    if terminal_oxygen == 1 {
        return Some(67);
    }
    if terminal_oxygen >= 2 {
        return Some(45);
    }
    let mut double_bonded_cn = false;
    for neighbor in neighbors {
        if neighbor.order != BondOrder::Double {
            continue;
        }
        let element = mol.atom(neighbor.neighbor).element;
        double_bonded_cn = matches!(element, Element::N | Element::C);
        if element == Element::C {
            for bond in bonds_of(mol, neighbor.neighbor) {
                if !double_bonded_cn {
                    break;
                }
                if bond.neighbor != idx {
                    double_bonded_cn = !(mol.atom(bond.neighbor).element == Element::N
                        && total_degree(mol, bond.neighbor) == 3);
                }
            }
        }
    }
    double_bonded_cn.then_some(54)
}

#[derive(Default)]
struct NitrogenConjugation {
    nco_or_ncs: bool,
    ncn_plus: bool,
    ngd_plus: bool,
    nnn_or_nnc: bool,
    neighbor_carbon: bool,
    neighbor_benzene_carbon: bool,
    double_bonded_to_carbon: Option<Element>,
    triple_bonded_to_carbon: Option<Element>,
    oxygen_bonded_to_carbon: usize,
    sulfur_bonded_to_carbon: usize,
}

fn conjugated_degree_three_n_type(
    mol: &Molecule,
    rings: &[Vec<AtomIdx>],
    neighbors: &[BondInfo],
    nso_family: &mut bool,
) -> Option<u8> {
    let mut state = NitrogenConjugation::default();
    for neighbor in neighbors {
        match mol.atom(neighbor.neighbor).element {
            Element::C => state.observe_carbon(mol, rings, neighbor.neighbor),
            Element::N => state.observe_nitrogen(mol, neighbor.neighbor),
            _ => {}
        }
    }
    if state.neighbor_carbon {
        *nso_family |= state.triple_bonded_to_carbon == Some(Element::N);
        if state.ncn_plus {
            return Some(55);
        }
        if state.ngd_plus {
            return Some(56);
        }
        if state.is_conjugated_amine(*nso_family) {
            return Some(40);
        }
    }
    (!*nso_family && (state.nco_or_ncs || state.nnn_or_nnc)).then_some(10)
}

impl NitrogenConjugation {
    fn observe_carbon(&mut self, mol: &Molecule, rings: &[Vec<AtomIdx>], carbon: AtomIdx) {
        self.neighbor_carbon = true;
        self.neighbor_benzene_carbon |=
            mol.atom(carbon).aromatic && atom_in_ring_of_size(rings, carbon, 6);
        let mut degree_two_n = 0;
        let mut degree_three_n = 0;
        let mut charged_n = 0;
        let mut aromatic_six_n = 0;
        self.oxygen_bonded_to_carbon = 0;
        self.sulfur_bonded_to_carbon = 0;

        for bond in bonds_of(mol, carbon) {
            let atom = bond.neighbor;
            let element = mol.atom(atom).element;
            self.nco_or_ncs |=
                bond.order == BondOrder::Double && matches!(element, Element::O | Element::S);
            let ring_count = rings.iter().filter(|ring| ring.contains(&atom)).count();
            if bond.order == BondOrder::Double
                || (bond.order == BondOrder::Aromatic
                    && (element == Element::C || (element == Element::N && ring_count == 1)))
            {
                self.double_bonded_to_carbon = Some(element);
            }
            if bond.order == BondOrder::Triple {
                self.triple_bonded_to_carbon = Some(element);
            }
            if element == Element::N && total_degree(mol, atom) == 3 {
                charged_n += usize::from(mol.atom(atom).charge == 1);
                aromatic_six_n +=
                    usize::from(atom_in_aromatic_ring_of_size(mol, rings, carbon, 6));
                let oxygen_on_n = bonds_of(mol, atom)
                    .iter()
                    .filter(|next| mol.atom(next.neighbor).element == Element::O)
                    .count();
                degree_three_n += usize::from(oxygen_on_n < 2);
            }
            degree_two_n += usize::from(
                element == Element::N
                    && total_degree(mol, atom) == 2
                    && matches!(bond.order, BondOrder::Double | BondOrder::Aromatic),
            );
            if mol.atom(atom).aromatic {
                self.oxygen_bonded_to_carbon += usize::from(element == Element::O);
                self.sulfur_bonded_to_carbon += usize::from(element == Element::S);
            }
        }
        if self.double_bonded_to_carbon == Some(Element::N) {
            self.ncn_plus |= degree_three_n == 2
                && degree_two_n == 0
                && charged_n > 0
                && aromatic_six_n == 0
                && total_degree(mol, carbon) < 4;
            self.ngd_plus |= degree_three_n == 3;
        }
    }

    fn observe_nitrogen(&mut self, mol: &Molecule, nitrogen: AtomIdx) {
        let mut counts = [0usize; 3];
        for bond in bonds_of(mol, nitrogen) {
            if bond.order != BondOrder::Double {
                continue;
            }
            let atom = bond.neighbor;
            match mol.atom(atom).element {
                Element::C => {
                    for next in bonds_of(mol, atom) {
                        if next.neighbor == nitrogen {
                            continue;
                        }
                        match mol.atom(next.neighbor).element {
                            Element::N => counts[0] += 1,
                            Element::O => counts[1] += 1,
                            Element::S => counts[2] += 1,
                            _ => {}
                        }
                    }
                    self.nnn_or_nnc |=
                        counts == [0, 0, 0] && !self.neighbor_benzene_carbon;
                }
                Element::N => self.nnn_or_nnc |= !self.neighbor_benzene_carbon,
                _ => {}
            }
        }
    }

    fn is_conjugated_amine(&self, nso_family: bool) -> bool {
        !self.nco_or_ncs
            && !nso_family
            && ((self.oxygen_bonded_to_carbon == 0
                && self.sulfur_bonded_to_carbon == 0
                && self.neighbor_benzene_carbon)
                || matches!(
                    self.double_bonded_to_carbon,
                    Some(Element::C | Element::N | Element::P)
                )
                || self.triple_bonded_to_carbon == Some(Element::C))
    }
}

fn degree_two_n_type(
    mol: &Molecule,
    idx: AtomIdx,
    neighbors: &[BondInfo],
    total_bond_order: u32,
    terminal_oxygen: usize,
    nso_family: bool,
) -> Option<u8> {
    if total_bond_order == 4 {
        return Some(if neighbors.iter().any(|bond| bond.order == BondOrder::Triple) {
            61
        } else {
            53
        });
    }
    if total_bond_order == 3 {
        let mut nitroso = false;
        let mut imine_or_azo = false;
        for neighbor in neighbors {
            if mol
                .bond_between(idx, neighbor.neighbor)
                .is_some_and(|(_, bond)| bond.order == BondOrder::Double)
            {
                let element = mol.atom(neighbor.neighbor).element;
                nitroso = element == Element::O && terminal_oxygen == 1;
                imine_or_azo = matches!(element, Element::C | Element::N);
            }
        }
        if nitroso && !imine_or_azo {
            return Some(46);
        }
        if imine_or_azo {
            return Some(9);
        }
    }
    if total_bond_order < 2 {
        return None;
    }
    let nso = neighbors.iter().any(|neighbor| {
        mol.atom(neighbor.neighbor).element == Element::S
            && bonds_of(mol, neighbor.neighbor)
                .iter()
                .filter(|bond| {
                    mol.atom(bond.neighbor).element == Element::O
                        && total_degree(mol, bond.neighbor) == 1
                })
                .count()
                == 1
    });
    if nso {
        Some(48)
    } else if !nso_family {
        Some(62)
    } else {
        None
    }
}

fn degree_one_n_type(mol: &Molecule, idx: AtomIdx, neighbors: &[BondInfo]) -> Option<u8> {
    let sp = neighbors.iter().any(|neighbor| {
        mol.bond_between(idx, neighbor.neighbor)
            .is_some_and(|(_, bond)| bond.order == BondOrder::Triple)
    });
    if sp {
        return Some(42);
    }
    let azide = neighbors.iter().any(|neighbor| {
        mol.atom(neighbor.neighbor).element == Element::N
            && total_degree(mol, neighbor.neighbor) == 2
            && bonds_of(mol, neighbor.neighbor).iter().any(|next| {
                let atom = next.neighbor;
                (mol.atom(atom).element == Element::N && total_degree(mol, atom) == 2)
                    || (mol.atom(atom).element == Element::C && total_degree(mol, atom) == 3)
            })
    });
    azide.then_some(47)
}

/// Faithful port of RDKit's aromatic-nitrogen cases (`AtomTyper.cpp`
/// `setMMFFHeavyAtomType`, `case 7:` under both the 5-ring and 6-ring
/// aromatic switches). Returns `None` if the atom is flagged aromatic but
/// isn't actually in a fully-aromatic-bonded 5- or 6-membered ring.
///
/// Not ported: the 5-ring N-oxide alpha/beta/other sub-distinction (N5AX
/// vs N5BX vs N5OX) -- RDKit itself collapses all three to numeric type 82,
/// so this does too, faithfully, not as a simplification of RDKit's own
/// behavior.
fn aromatic_n_type(mol: &Molecule, rings: &[Vec<AtomIdx>], idx: AtomIdx) -> Option<u8> {
    if atom_in_aromatic_ring_of_size(mol, rings, idx, 5) {
        if is_atom_n_oxide(mol, idx) {
            return Some(82); // N5AX/N5BX/N5OX, collapsed like RDKit's own code
        }
        let het = find_alpha_beta_heteroatoms(mol, rings, idx);
        if het.alpha.is_empty() && het.beta.is_empty() {
            if total_degree(mol, idx) == 3 {
                return Some(39); // NPYL: pyrrole-type N with pi lone pair
            }
            return Some(76); // N5M: anionic 5-ring aromatic N
        }
        if total_degree(mol, idx) == 3 && het.alpha.len() != het.beta.len() {
            return Some(81); // NIM+/N5A+/N5B+/N5+: positively charged 5-ring N
        }
        if !het.alpha.is_empty() && (het.beta.is_empty() || het.is_alpha_os) {
            return Some(65); // N5A: alpha to N/O/S
        }
        if !het.beta.is_empty() && (het.alpha.is_empty() || het.is_beta_os) {
            return Some(66); // N5B: beta to N/O/S
        }
        if !het.alpha.is_empty() && !het.beta.is_empty() {
            return Some(79); // N5: general 5-ring aromatic N
        }
    }

    if atom_in_aromatic_ring_of_size(mol, rings, idx, 6) {
        if is_atom_n_oxide(mol, idx) {
            return Some(69); // NPOX: pyridinium N-oxide
        }
        if total_degree(mol, idx) == 3 {
            return Some(58); // NPD+: protonated pyridinium N
        }
        return Some(38); // NPYD: neutral pyridine-type N
    }

    None
}

// ── O type assignment ────────────────────────────────────────────────────────

/// A genuinely terminal oxygen in RDKit's sense (`AtomTyper.cpp`'s
/// `atom->getDegree() <= 1` gate at the pinned commit's line 1554, reached
/// only after the 3- and (total-degree-)2-neighbor branches have already
/// been ruled out): at most one heavy neighbor *and* no implicit H. This
/// excludes -OH (1 heavy neighbor + 1 implicit H, `getTotalDegree()==2`)
/// and bridging ether/ester -O- (2 heavy neighbors) alike, leaving only
/// =O / -O⁻-style oxygens whose valence is already fully satisfied by
/// their one explicit bond.
fn is_terminal_o(mol: &Molecule, idx: AtomIdx) -> bool {
    mol.atom(idx).element == Element::O
        && bonds_of(mol, idx).len() <= 1
        && implicit_hcount(mol, idx) == 0
}

/// RDKit's `nObondedToCorNorS`, O contribution: count of terminal oxygens
/// bonded to `central` (`AtomTyper.cpp` lines 1597-1600). Note this counts
/// `central`'s terminal-O neighbors generally, so when called with
/// `central` = the ipso oxygen's own neighbor, it naturally includes the
/// ipso oxygen itself alongside any sibling terminal oxygens -- exactly as
/// RDKit's loop does (it never excludes the atom it was reached from).
fn count_terminal_o_neighbors(mol: &Molecule, central: AtomIdx) -> usize {
    bonds_of(mol, central)
        .iter()
        .filter(|b| is_terminal_o(mol, b.neighbor))
        .count()
}

/// RDKit's `nNbondedToCorNorS`: count of degree-2 (`getTotalDegree()==2`)
/// nitrogens bonded to `central` (`AtomTyper.cpp` lines 1593-1596).
fn count_deg2_n_neighbors(mol: &Molecule, central: AtomIdx) -> usize {
    bonds_of(mol, central)
        .iter()
        .filter(|b| {
            mol.atom(b.neighbor).element == Element::N && total_degree(mol, b.neighbor) == 2
        })
        .count()
}

/// RDKit's `nSbondedToCorNorS`: count of terminal (`getTotalDegree()==1`)
/// sulfurs bonded to `central` (`AtomTyper.cpp` lines 1601-1604).
fn count_terminal_s_neighbors(mol: &Molecule, central: AtomIdx) -> usize {
    bonds_of(mol, central)
        .iter()
        .filter(|b| {
            mol.atom(b.neighbor).element == Element::S && total_degree(mol, b.neighbor) == 1
        })
        .count()
}

/// Source-grounded deterministic port of RDKit's `case 8:` degree-≤1
/// ("terminal") oxygen branch (`AtomTyper.cpp` lines 1554-1737 at the
/// pinned commit -- issue #227 Priority 1A-3). Distinguishes OM (35, oxide
/// oxygen on an unremarkable sp3 C/N or bonded to literal H), O2CM (32, the
/// carboxylate / nitro-nitrate / genuine N-oxide / thiosulfinate /
/// sulfate-sulfonate-sulfonamide-sulfone / phosphate-phosphonate-
/// phosphine-oxide / perchlorate terminal-oxygen union), and O=C (7,
/// generic carbonyl/nitroso/sulfoxide oxygen) for a genuinely terminal
/// oxygen. Returns `None` if `idx` is not terminal or if none of RDKit's
/// conditions fire (falls back to `assign_o_type`'s pre-existing generic
/// handling below).
///
/// This is not a "central element implies type" shortcut: routing depends
/// on this bond's order and on the formal degree/valence of the *central*
/// atom and how many *other* terminal O/S atoms (or degree-2 N atoms) share
/// it, exactly as RDKit computes it -- e.g. a carboxylic acid's neutral
/// C=O is *not* O2CM (its carbon has only one terminal oxygen), only a
/// carboxylate *anion*'s two terminal oxygens are (both share a carbon
/// with two terminal oxygens), matching real MMFF94/RDKit output.
///
/// One RDKit quirk ported faithfully, not smoothed over: for a phosphorus
/// or chlorine central atom, `isPhosphateOrPerchlorateO` really is
/// unconditional in RDKit's own source (`atomicNum==15 || atomicNum==17`,
/// no bond-order or valence check at all) -- this is RDKit's own
/// documented behavior (its type-32 symbol list explicitly names "OP,
/// Oxygen in phosphine oxide"), not a chematic simplification, so a
/// neutral phosphine oxide's O also becomes O2CM here, same as RDKit.
fn classify_terminal_o(mol: &Molecule, idx: AtomIdx) -> Option<u8> {
    let central_bond = bonds_of(mol, idx).into_iter().next()?;
    let central = central_bond.neighbor;
    let bond_order = central_bond.order;
    let central_elem = mol.atom(central).element;

    let n_o = count_terminal_o_neighbors(mol, central);
    let n_n = count_deg2_n_neighbors(mol, central);
    let n_s = count_terminal_s_neighbors(mol, central);

    let mut is_oxide_on_c_or_n = false;
    let mut is_o2cm = false;
    let mut is_carbonyl_like = false;

    match central_elem {
        Element::C => {
            if n_o == 2 {
                is_o2cm = true; // isCarboxylateO
            }
            if bond_order == BondOrder::Double {
                is_carbonyl_like = true; // isCarbonylO
            } else if bond_order == BondOrder::Single && n_o == 1 {
                is_oxide_on_c_or_n = true; // isOxideOBondedToC
            }
        }
        Element::N => {
            if n_o >= 2 {
                is_o2cm = true; // isNitroO (nitro/nitrate)
            }
            if bond_order == BondOrder::Double {
                is_carbonyl_like = true; // isNitrosoO
            } else if bond_order == BondOrder::Single && n_o == 1 {
                // RDKit distinguishes a genuine N-oxide (isNOxideO, ipso
                // -> O2CM) from a plain oxide oxygen on an otherwise-normal
                // trivalent nitrogen (isOxideOBondedToN, ipso -> OM) by
                // whether the *central* N's real bond-order-sum valence is
                // 4 vs 2-or-3. Summing bond orders directly is unsafe here:
                // `central` may be a ring atom whose bonds are still
                // `BondOrder::Aromatic` in this post-`compute_mmff94_aromatic_view`
                // molecule (order_int()==1, undercounting a ring bond that's
                // really a Kekule double bond -- e.g. pyridine N-oxide's
                // ring N would wrongly read total_valence==3, not 4). Formal
                // charge is a bond-order-representation-independent proxy
                // that holds for every standard, RDKit-sanitizable N-oxide
                // depiction (charge-separated N+/O-, the only form RDKit
                // itself produces/accepts): a genuinely oxidized nitrogen
                // carries a +1 formal charge; a merely-oxide-substituted
                // normal trivalent nitrogen does not. A non-charge-separated
                // "neutral pentavalent N=O" spelling, if it parses at all,
                // is out of scope for this port (RDKit's own sanitizer does
                // not produce or expect that form).
                if mol.atom(central).charge > 0 {
                    is_o2cm = true; // isNOxideO: genuine N-oxide
                } else {
                    is_oxide_on_c_or_n = true; // isOxideOBondedToN
                }
            }
        }
        Element::S => {
            if n_s == 1 {
                is_o2cm = true; // isThioSulfinateO
            }
            if bond_order == BondOrder::Single
                || (bond_order == BondOrder::Double && (n_o + n_n) > 1)
            {
                is_o2cm = true; // isSulfateO: sulfate/sulfonate/sulfonamide/sulfone
            } else if bond_order == BondOrder::Double && (n_o + n_n) == 1 {
                is_carbonyl_like = true; // isSulfoxideO
            }
        }
        Element::P | Element::CL => {
            is_o2cm = true; // isPhosphateOrPerchlorateO -- unconditional, see doc above
        }
        Element::H => {
            is_oxide_on_c_or_n = true; // isOxideOBondedToH
        }
        _ => return None,
    }

    if is_oxide_on_c_or_n {
        return Some(35); // OM
    }
    if is_o2cm {
        return Some(32); // O2CM
    }
    if is_carbonyl_like {
        return Some(7); // O=C / O=N / O=S generic
    }
    None
}

fn assign_o_type(
    mol: &Molecule,
    rings: &[Vec<AtomIdx>],
    idx: AtomIdx,
) -> Result<u8, NumericTypeError> {
    // Aromatic 5-ring oxygen (furan-type) -- RDKit's aromatic switch has no
    // 6-ring case for oxygen (no common neutral 6-ring aromatic O in the
    // chemistry this typer covers), so only the 5-ring case is ported.
    if mol.atom(idx).aromatic && atom_in_aromatic_ring_of_size(mol, rings, idx, 5) {
        return Ok(59); // OFUR
    }

    // Terminal oxygen (=O, -O⁻, ...): issue #227 Priority 1A-3's OM(35) vs
    // O2CM(32) vs O=C(7) disambiguation. Falls through to the generic
    // handling below for non-terminal O (ether/ester/alcohol) or the rare
    // terminal O whose central atom isn't C/N/S/P/Cl/H.
    if is_terminal_o(mol, idx)
        && let Some(t) = classify_terminal_o(mol, idx)
    {
        return Ok(t);
    }

    // RDKit `case 8` for two- and three-connected oxygen.
    match total_degree(mol, idx) {
        3 => return Ok(49), // O+: oxonium
        2 => {
            // (An MMFF-aromatic ring's bonds count 1 here; a two-connected
            // O+ has valence three either way.)
            if total_valence(mol, idx) == 3 || mol.atom(idx).charge == 1 {
                return Ok(51); // O=+: oxenium (pyrylium)
            }
            let h = bonds_of(mol, idx)
                .iter()
                .filter(|b| mol.atom(b.neighbor).element == Element::H)
                .count()
                + usize::from(implicit_hcount(mol, idx));
            if h == 2 {
                return Ok(70); // OH2: water
            }
            return Ok(6); // OR
        }
        _ => {}
    }

    // Double bond to C or N → carbonyl/similar oxygen (type 7)
    if count_bond_order(mol, idx, BondOrder::Double) > 0 {
        return Ok(7); // O=C
    }

    // Anionic O (formal charge -1) not resolved above: phenoxide and
    // similar oxide oxygens on an unremarkable carbon/nitrogen.
    // Registry-verified: type 35 is OM (OXIDE OXYGEN ON SP3 C); type 34
    // is NR+, a nitrogen-only type -- the previous `34` here was the
    // same class of silent element collision fixed in `assign_n_type`
    // above (that function's protonated-N branch had the mirror-image
    // bug: it returned 32, O2CM's id, instead of 34).
    if mol.atom(idx).charge < 0 {
        return Ok(35); // OM
    }

    // Single-bond O (ether, alcohol, ester, amide O)
    Ok(6) // OR
}

// ── S type assignment ────────────────────────────────────────────────────────

/// Source-grounded deterministic port of RDKit's `case 16:` (sulfur)
/// degree-based branch (`AtomTyper.cpp` lines ~1815-1917 at the pinned
/// commit, issue #227). The previous version only counted *explicit double
/// bonds* to oxygen, which missed every charge-separated sulfoxide/sulfone
/// spelling (e.g. `[S+]([O-])`, MMFF94's only valid form for a charged
/// sulfinyl/sulfonyl sulfur) -- this port instead follows RDKit's own
/// degree/terminal-neighbor-counting rule, which is bond-order-separation
/// agnostic by construction (a terminal O counts whether it's reached via a
/// double bond or a charge-separated single bond).
fn assign_s_type(
    mol: &Molecule,
    rings: &[Vec<AtomIdx>],
    idx: AtomIdx,
) -> Result<u8, NumericTypeError> {
    let atom = mol.atom(idx);
    // S5 is the sulfur of a five-membered aromatic ring; an aromatic sulfur
    // in a six-membered ring (methylene blue's `[S+]`) takes the degree
    // rules below, as in RDKit.
    if atom.aromatic && rings.iter().any(|r| r.len() == 5 && r.contains(&idx)) {
        return Ok(44); // S5 aromatic sulfur (thiophene)
    }

    let degree = total_degree(mol, idx);
    let nbrs = bonds_of(mol, idx);

    // 3- or 4-connected sulfur (sulfoxide/sulfone/thiosulfinate family).
    if degree == 3 || degree == 4 {
        return Ok(multicoordinate_s_type(mol, &nbrs, degree));
    }

    // 2-connected sulfur.
    if degree == 2 {
        let o_double_bonded = nbrs
            .iter()
            .any(|b| b.order == BondOrder::Double && mol.atom(b.neighbor).element == Element::O);
        if o_double_bonded {
            return Ok(74); // =S=O: sulfinyl sulfur (e.g. C=S=O)
        }
        return Ok(15); // S: thiol, sulfide, or disulfide
    }

    // 1-connected (terminal) sulfur.
    if degree == 1
        && let Some(nb) = nbrs.first()
    {
        return Ok(terminal_s_type(mol, nb));
    }

    Ok(15) // conservative fallback (degree 0, e.g. an isolated ion)
}

fn multicoordinate_s_type(mol: &Molecule, neighbors: &[BondInfo], degree: usize) -> u8 {
    let terminal_oxygen_or_nitrogen = neighbors
        .iter()
        .filter(|bond| {
            let atom = mol.atom(bond.neighbor);
            (atom.element == Element::O && bonds_of(mol, bond.neighbor).len() == 1)
                || (atom.element == Element::N && total_degree(mol, bond.neighbor) == 2)
        })
        .count();
    let terminal_sulfur = neighbors
        .iter()
        .filter(|bond| {
            mol.atom(bond.neighbor).element == Element::S
                && bonds_of(mol, bond.neighbor).len() == 1
        })
        .count();
    let double_bonded_carbon = neighbors.iter().any(|bond| {
        bond.order == BondOrder::Double && mol.atom(bond.neighbor).element == Element::C
    });
    if (degree == 3 && terminal_oxygen_or_nitrogen == 2 && double_bonded_carbon) || degree == 4 {
        18
    } else if (terminal_oxygen_or_nitrogen > 0 && terminal_sulfur > 0)
        || (terminal_oxygen_or_nitrogen == 2 && !double_bonded_carbon)
    {
        73
    } else {
        17
    }
}

fn terminal_s_type(mol: &Molecule, neighbor: &BondInfo) -> u8 {
    let terminal_sulfur = bonds_of(mol, neighbor.neighbor)
        .iter()
        .filter(|bond| {
            mol.atom(bond.neighbor).element == Element::S
                && bonds_of(mol, bond.neighbor).len() == 1
        })
        .count();
    let double_bonded_carbon = neighbor.order == BondOrder::Double
        && mol.atom(neighbor.neighbor).element == Element::C;
    if double_bonded_carbon && terminal_sulfur != 2 {
        16
    } else {
        72
    }
}

// ── P type assignment ────────────────────────────────────────────────────────

fn assign_p_type(mol: &Molecule, idx: AtomIdx) -> Result<u8, NumericTypeError> {
    // RDKit `setMMFFHeavyAtomType` `case 15`: by total degree only. (The
    // previous "P=O -> 25" rule typed a thiophosphate P(=S)(O)(O)O as 26.)
    Ok(match total_degree(mol, idx) {
        4 => 25, // PO4 / PO3 / PO2 / PO / PTET
        2 => 75, // -P=C
        _ => 26, // P: tricoordinate phosphorus
    })
}

// ── H type assignment ────────────────────────────────────────────────────────

/// MMFF94 type of an explicit hydrogen, derived from its parent atom's
/// already-assigned MMFF type the way RDKit's typer (and MMFF's `MMFFHDEF.PAR`)
/// does: the parent's symbolic type decides the hydrogen type, with the
/// oxygen-bound hydrogen of a generic `OR` oxygen (type 6) refined by what
/// the oxygen is attached to.
fn assign_h_type(mol: &Molecule, idx: AtomIdx, types: &[u8]) -> Result<u8, NumericTypeError> {
    let nbrs = bonds_of(mol, idx);
    if nbrs.is_empty() {
        return Ok(5); // H_C fallback
    }
    let parent = nbrs[0].neighbor;
    let parent_type = types[parent.0 as usize];
    let nbr_atom = mol.atom(parent);

    Ok(match nbr_atom.element {
        Element::C | Element::SI => 5, // HC / HSI
        Element::O => assign_oxygen_bound_h_type(mol, parent, parent_type),
        Element::S => 71, // HS  H on sulfur
        Element::P => 71, // HP  H on phosphorus
        Element::N => match parent_type {
            8 | 39 | 62 => 23,                 // HNR / HPYL / H on NR-
            9 => 27,                           // HN=C / HN=N
            10 | 40 | 43 => 28,                // HNCO, HNCS, HNCC, HNCN, HNSO
            34 | 54 | 55 | 56 | 58 | 81 => 36, // HNR+, HN=C+, HNCN+, HGD+, HPD+, HIM+
            _ => legacy_nitrogen_bound_h_type(mol, parent),
        },
        _ => 5,
    })
}

/// Local-environment N-H typing used only for parent types without an
/// `MMFFHDEF` row above (rare N-oxide and nitroso-like parents).
fn legacy_nitrogen_bound_h_type(mol: &Molecule, n_idx: AtomIdx) -> u8 {
    let n_atom = mol.atom(n_idx);
    if n_atom.aromatic {
        return 23; // treat as HNR for aromatic NH
    }
    let n_is_amide = bonds_of(mol, n_idx).iter().any(|b| {
        b.order == BondOrder::Single
            && mol.atom(b.neighbor).element == Element::C
            && bonds_of(mol, b.neighbor).iter().any(|bb| {
                bb.order == BondOrder::Double && mol.atom(bb.neighbor).element == Element::O
            })
    });
    if n_is_amide {
        28 // HNCO H on amide N
    } else if count_bond_order(mol, n_idx, BondOrder::Double) > 0 {
        27 // HN=C H on imine N
    } else {
        23 // HNR  H on amine N
    }
}

/// Select the MMFF94 numeric type for an explicit hydrogen bonded to oxygen.
///
/// Matches RDKit 2026.03.6 (the pinned oracle; probe table in
/// `mmff_oxygen_h_types_follow_rdkit`):
///
/// * water H is HOH (31), hydroxide H is HOM (21);
/// * oxonium parents give HO+ (50) and HO=+ (52);
/// * on a generic `OR` oxygen the H is HOS (33) when the oxygen is bonded to
///   sulfur, HOP (24) when bonded to phosphorus, HOCO (24) when bonded to a
///   carbonyl carbon, HOCC (29) when bonded to an enol/imidic sp2 carbon
///   (C=C, C=N) or an aromatic carbon, and HOR (21) otherwise — including
///   N-OH (hydroxylamines, oximes, hydroxamic acids) and thioacid O-H.
fn assign_oxygen_bound_h_type(mol: &Molecule, oxygen_idx: AtomIdx, oxygen_type: u8) -> u8 {
    match oxygen_type {
        49 => return 50, // HO+
        51 => return 52, // HO=+
        70 => return 31, // HOH
        _ => {}
    }
    let heavy_neighbors: Vec<_> = bonds_of(mol, oxygen_idx)
        .into_iter()
        .filter(|bond| mol.atom(bond.neighbor).element != Element::H)
        .collect();

    if heavy_neighbors.is_empty() {
        let hydrogens = bonds_of(mol, oxygen_idx).len();
        return if hydrogens >= 2 { 31 } else { 21 }; // HOH water / HOM hydroxide
    }

    let heavy_neighbor = heavy_neighbors[0].neighbor;
    let heavy_atom = mol.atom(heavy_neighbor);
    match heavy_atom.element {
        Element::S => 33, // HOS
        Element::P => 24, // HOP
        Element::C => {
            let carbon_bonds = bonds_of(mol, heavy_neighbor);
            let double_to = |element: Element| {
                carbon_bonds.iter().any(|bond| {
                    bond.order == BondOrder::Double && mol.atom(bond.neighbor).element == element
                })
            };
            if double_to(Element::O) {
                24 // HOCO, acid hydroxyl
            } else if heavy_atom.aromatic || double_to(Element::C) || double_to(Element::N) {
                29 // HOCC, enol / phenol / imidic acid hydroxyl
            } else {
                21 // HOR
            }
        }
        _ => 21, // HOR, including N-OH
    }
}
