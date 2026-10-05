//! RDKit 2026.03.6 reaction stereochemistry, for the pinned RDKit profile
//! (issue #734).
//!
//! RDKit does not use a reactant template's `@`/`@@` when matching, and it
//! sets a product atom's chiral tag from tags and *bond orders* rather than
//! from geometry (`ReactionUtils.cpp::updateProductsStereochem`,
//! `ReactionRunner.cpp::checkProductChirality` and
//! `checkAndCorrectChiralityOfMatchingAtomsInProduct`):
//!
//! - Each mapped product-template atom gets an inversion flag by comparing
//!   the raw tags of the two template atoms: 1 invert, 2 retain, 3 remove,
//!   4 set from the product template (also every unmapped atom).
//! - The product atom first takes the reactant atom's raw tag (flag 1 or 2),
//!   referred to the product atom's bond order — not the reactant's.
//! - When the product atom's neighbours can be traced back to the reactant
//!   atom's bonds (degree 3+, degrees within one, at most one unknown
//!   neighbour), the tag is recomputed from the permutation between the two
//!   bond orders, so the reactant's configuration carries over (inverted for
//!   flag 1). Otherwise the raw copy stands.
//!
//! The raw copy is why RDKit gives D-alanine for `C[C@H](N)C(=O)O` under the
//! identity template `[N:1][C@@H:2](C)C(=O)O>>[N:1][C@@H:2](C)C(=O)O`. To
//! reproduce it this module rebuilds RDKit's bond order: for templates from
//! their text (chain bonds in creation order, then ring closures by ring
//! number), for reactants from the chematic SMILES parser's bond creation
//! order. The one thing a parsed molecule does not record is the ring-closure
//! *number*, which orders an atom's ring bonds in RDKit; when an atom has two
//! or more ring closures every order is tried, and a product whose stereo
//! depends on the choice is reported as unsupported instead of guessed.

use chematic_core::{AtomIdx, Chirality, Molecule, STEREO_H_SENTINEL};
use rustc_hash::{FxHashMap, FxHashSet};

/// Largest number of bond-order alternatives tried for one product atom.
const MAX_ORDER_OPTIONS: usize = 720;

/// One template atom as RDKit's SMARTS parser builds it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct TemplateAtom {
    pub map: Option<u16>,
    /// Chiral tag relative to `order` (and an implicit H last): `Some(true)`
    /// is `CHI_TETRAHEDRAL_CW` (`@@`), `Some(false)` `CHI_TETRAHEDRAL_CCW`.
    pub tag: Option<bool>,
    /// Neighbouring template atoms in RDKit bond order.
    pub order: Vec<usize>,
}

/// Stereo facts for one prepared reaction.
#[derive(Debug, Clone)]
pub(crate) struct ReactionStereo {
    pub reactants: Vec<Vec<TemplateAtom>>,
    pub products: Vec<Vec<TemplateAtom>>,
    /// Per product template, per atom: RDKit's `molInversionFlag` (1–4).
    pub flags: Vec<Vec<Option<u8>>>,
    /// Map-number pairs bonded in a reactant template (unmapped = 0).
    pub reactant_bond_maps: FxHashSet<(u16, u16)>,
}

/// Why a product's stereo cannot be reproduced.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum StereoUnsupported {
    /// The result depends on a bond order the reactant does not record.
    BondOrder,
}

impl ReactionStereo {
    /// Build from the reactant and product sides of one (normalized) SMIRKS.
    /// `None` when a template uses syntax this reading does not follow
    /// (component grouping, non-tetrahedral chirality classes, `@?`).
    pub fn new(reactant_side: &str, product_side: &str) -> Option<Self> {
        let reactants = split_templates(reactant_side)
            .into_iter()
            .map(tokenize_template)
            .collect::<Option<Vec<_>>>()?;
        let products = split_templates(product_side)
            .into_iter()
            .map(tokenize_template)
            .collect::<Option<Vec<_>>>()?;
        let flags = inversion_flags(&reactants, &products);
        let mut reactant_bond_maps = FxHashSet::default();
        for template in &reactants {
            for atom in template {
                for &n in &atom.order {
                    let a = atom.map.unwrap_or(0);
                    let b = template[n].map.unwrap_or(0);
                    reactant_bond_maps.insert((a.min(b), a.max(b)));
                }
            }
        }
        Some(Self {
            reactants,
            products,
            flags,
            reactant_bond_maps,
        })
    }
}

/// Top-level `.`-separated templates of one reaction side.
fn split_templates(side: &str) -> Vec<&str> {
    let mut parts = Vec::new();
    let mut depth = 0usize;
    let mut start = 0;
    for (i, b) in side.bytes().enumerate() {
        match b {
            b'[' => depth += 1,
            b']' => depth = depth.saturating_sub(1),
            b'.' if depth == 0 => {
                parts.push(&side[start..i]);
                start = i + 1;
            }
            _ => {}
        }
    }
    parts.push(&side[start..]);
    parts.into_iter().filter(|p| !p.is_empty()).collect()
}

/// One entry of an atom's SMILES-text neighbour order.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Entry {
    Atom(usize),
    H,
    Pending,
}

struct Building {
    map: Option<u16>,
    /// 0 none, 1 `@`, 2 `@@`.
    at: u8,
    h_adjust: bool,
    text: Vec<Entry>,
}

/// Parse one template component the way RDKit's SMARTS parser orders bonds
/// and chiral tags.
pub(crate) fn tokenize_template(text: &str) -> Option<Vec<TemplateAtom>> {
    let b = text.as_bytes();
    let mut atoms: Vec<Building> = Vec::new();
    let mut chain: Vec<(usize, usize)> = Vec::new();
    let mut rings: Vec<((u32, u32), usize, usize)> = Vec::new();
    let mut open: FxHashMap<u32, (usize, usize)> = FxHashMap::default();
    let mut occurrences: FxHashMap<u32, u32> = FxHashMap::default();
    let mut prev: Option<usize> = None;
    let mut stack: Vec<Option<usize>> = Vec::new();
    let mut i = 0;
    let add_atom = |atoms: &mut Vec<Building>,
                    prev: &mut Option<usize>,
                    chain: &mut Vec<(usize, usize)>,
                    mut atom: Building| {
        let x = atoms.len();
        if let Some(p) = *prev {
            chain.push((p, x));
            atom.text.push(Entry::Atom(p));
            atoms[p].text.push(Entry::Atom(x));
        }
        if atom.h_adjust {
            atom.text.push(Entry::H);
        }
        atoms.push(atom);
        *prev = Some(x);
    };
    while i < b.len() {
        match b[i] {
            b'[' => {
                let mut depth = 0usize;
                let mut close = None;
                for (j, &c) in b.iter().enumerate().skip(i) {
                    match c {
                        b'[' => depth += 1,
                        b']' => {
                            depth -= 1;
                            if depth == 0 {
                                close = Some(j);
                                break;
                            }
                        }
                        _ => {}
                    }
                }
                let close = close?;
                let atom = bracket_atom(&text[i + 1..close])?;
                add_atom(&mut atoms, &mut prev, &mut chain, atom);
                i = close + 1;
            }
            b'(' => {
                stack.push(prev);
                i += 1;
            }
            b')' => {
                prev = stack.pop()?;
                i += 1;
            }
            b'.' => {
                // Component grouping inside one template is not followed.
                return None;
            }
            c @ (b'0'..=b'9' | b'%') => {
                let (num, len) = if c == b'%' {
                    if b.get(i + 1) == Some(&b'(') {
                        let end = text[i + 2..].find(')')? + i + 2;
                        (text[i + 2..end].parse::<u32>().ok()?, end + 1 - i)
                    } else {
                        (text.get(i + 1..i + 3)?.parse::<u32>().ok()?, 3)
                    }
                } else {
                    (u32::from(c - b'0'), 1)
                };
                let x = prev?;
                if let Some((y, slot)) = open.remove(&num) {
                    let occurrence = occurrences.entry(num).or_insert(0);
                    rings.push(((num, *occurrence), y, x));
                    *occurrence += 1;
                    atoms[x].text.push(Entry::Atom(y));
                    atoms[y].text[slot] = Entry::Atom(x);
                } else {
                    open.insert(num, (x, atoms[x].text.len()));
                    atoms[x].text.push(Entry::Pending);
                }
                i += len;
            }
            b'C' if b.get(i + 1) == Some(&b'l') => {
                add_atom(&mut atoms, &mut prev, &mut chain, plain_atom());
                i += 2;
            }
            b'B' if b.get(i + 1) == Some(&b'r') => {
                add_atom(&mut atoms, &mut prev, &mut chain, plain_atom());
                i += 2;
            }
            b'B' | b'C' | b'N' | b'O' | b'P' | b'S' | b'F' | b'I' | b'b' | b'c' | b'n' | b'o'
            | b'p' | b's' | b'a' | b'A' | b'*' => {
                add_atom(&mut atoms, &mut prev, &mut chain, plain_atom());
                i += 1;
            }
            b'-' | b'=' | b'#' | b':' | b'~' | b'/' | b'\\' | b'$' | b'@' | b'!' | b';' | b','
            | b'&' | b'<' | b'>' => i += 1,
            _ => return None,
        }
    }
    if !open.is_empty() {
        return None;
    }
    // RDKit bond indices: chain bonds as created, then ring closures by ring
    // number and occurrence (`SmilesParseOps::CloseMolRings`).
    rings.sort_by_key(|&(key, _, _)| key);
    let mut neighbours: Vec<Vec<usize>> = vec![Vec::new(); atoms.len()];
    for (x, y) in chain
        .into_iter()
        .chain(rings.into_iter().map(|(_, y, x)| (y, x)))
    {
        neighbours[x].push(y);
        neighbours[y].push(x);
    }
    let mut out = Vec::with_capacity(atoms.len());
    for (atom, order) in atoms.into_iter().zip(neighbours) {
        let tag = match atom.at {
            0 => None,
            at => {
                let code = |e: &Entry| match *e {
                    Entry::Atom(n) => n as u32,
                    _ => STEREO_H_SENTINEL,
                };
                let text: Vec<u32> = atom.text.iter().map(code).collect();
                let mut rd: Vec<u32> = order.iter().map(|&n| n as u32).collect();
                if atom.h_adjust {
                    rd.push(STEREO_H_SENTINEL);
                }
                let even = permutation_parity(&text, &rd)?;
                Some((at == 2) == even)
            }
        };
        out.push(TemplateAtom {
            map: atom.map,
            tag,
            order,
        });
    }
    Some(out)
}

fn plain_atom() -> Building {
    Building {
        map: None,
        at: 0,
        h_adjust: false,
        text: Vec::new(),
    }
}

/// Map number, chirality and whether RDKit reads the bracket's `H` as the
/// implicit hydrogen of a stereocentre (only in the simple SMILES-like
/// spelling, `[C@@H:2]`; not in `[C@@;H1:2]`).
fn bracket_atom(inner: &str) -> Option<Building> {
    // Strip recursive SMARTS before looking for `@` and the map number.
    let mut flat = String::with_capacity(inner.len());
    let mut depth = 0usize;
    let bytes = inner.as_bytes();
    let mut k = 0;
    while k < bytes.len() {
        if depth == 0 && bytes[k] == b'$' && bytes.get(k + 1) == Some(&b'(') {
            depth = 1;
            k += 2;
            continue;
        }
        if depth > 0 {
            match bytes[k] {
                b'(' => depth += 1,
                b')' => depth -= 1,
                _ => {}
            }
            k += 1;
            continue;
        }
        flat.push(bytes[k] as char);
        k += 1;
    }
    let (body, map) = match flat.rfind(':') {
        Some(pos)
            if pos + 1 < flat.len() && flat[pos + 1..].bytes().all(|c| c.is_ascii_digit()) =>
        {
            (&flat[..pos], flat[pos + 1..].parse::<u16>().ok())
        }
        _ => (flat.as_str(), None),
    };
    let fb = body.as_bytes();
    let first_at = fb.iter().position(|&c| c == b'@');
    let at = match first_at {
        None => 0,
        Some(p) => {
            let n = if fb.get(p + 1) == Some(&b'@') { 2 } else { 1 };
            let after = fb.get(p + n as usize).copied();
            if fb[p + n as usize..].contains(&b'@')
                || matches!(after, Some(b'?' | b'T' | b'S' | b'O' | b'A'))
            {
                return None;
            }
            n
        }
    };
    let h_adjust = at > 0 && simple_chiral_h(body);
    Some(Building {
        map,
        at,
        h_adjust,
        text: Vec::new(),
    })
}

/// `isotope? symbol @{1,2} H<digit>? charge?` with an `H`.
fn simple_chiral_h(body: &str) -> bool {
    let b = body.as_bytes();
    let mut i = 0;
    while i < b.len() && b[i].is_ascii_digit() {
        i += 1;
    }
    if i < b.len() && b[i] == b'#' {
        i += 1;
        let start = i;
        while i < b.len() && b[i].is_ascii_digit() {
            i += 1;
        }
        if i == start {
            return false;
        }
    } else if i < b.len() && (b[i].is_ascii_alphabetic() || b[i] == b'*') {
        i += 1;
        if i < b.len() && b[i].is_ascii_lowercase() && b[i] != b'h' {
            i += 1;
        }
    } else {
        return false;
    }
    let at_start = i;
    while i < b.len() && b[i] == b'@' {
        i += 1;
    }
    if i == at_start || i < b.len() && b[i] != b'H' {
        return false;
    }
    if i == b.len() {
        return false;
    }
    i += 1; // H
    if i < b.len() && b[i].is_ascii_digit() {
        i += 1;
    }
    while i < b.len() && matches!(b[i], b'+' | b'-') {
        i += 1;
    }
    while i < b.len() && b[i].is_ascii_digit() {
        i += 1;
    }
    i == b.len()
}

/// Parity of the permutation taking `from` to `to`: `Some(true)` even.
pub(crate) fn permutation_parity(from: &[u32], to: &[u32]) -> Option<bool> {
    if from.len() != to.len() {
        return None;
    }
    let mut perm = Vec::with_capacity(to.len());
    for t in to {
        perm.push(from.iter().position(|f| f == t)?);
    }
    let mut seen = vec![false; perm.len()];
    for &p in &perm {
        if std::mem::replace(&mut seen[p], true) {
            return None;
        }
    }
    let mut inversions = 0usize;
    for i in 0..perm.len() {
        for j in i + 1..perm.len() {
            if perm[i] > perm[j] {
                inversions += 1;
            }
        }
    }
    Some(inversions.is_multiple_of(2))
}

/// RDKit's `countSwapsToInterconvert(ref, probe)`; `None` where RDKit would
/// fail its invariant (an element of `reference` missing from `probe`).
fn count_swaps<T: PartialEq + Copy>(reference: &[T], probe: &[T]) -> Option<usize> {
    if reference.len() != probe.len() {
        return None;
    }
    let mut probe = probe.to_vec();
    let mut swaps = 0;
    for i in 0..reference.len() {
        if probe[i] != reference[i] {
            let j = (i..probe.len()).find(|&j| probe[j] == reference[i])?;
            probe.swap(i, j);
            swaps += 1;
        }
    }
    Some(swaps)
}

/// `ReactionUtils.cpp::updateProductsStereochem`.
fn inversion_flags(
    reactants: &[Vec<TemplateAtom>],
    products: &[Vec<TemplateAtom>],
) -> Vec<Vec<Option<u8>>> {
    // The last reactant-template atom carrying each map number.
    let mut by_map: FxHashMap<u16, (usize, usize)> = FxHashMap::default();
    for (r, template) in reactants.iter().enumerate() {
        for (i, atom) in template.iter().enumerate() {
            if let Some(m) = atom.map {
                by_map.insert(m, (r, i));
            }
        }
    }
    products
        .iter()
        .map(|template| {
            template
                .iter()
                .map(|p_atom| {
                    let Some(m) = p_atom.map else {
                        return Some(4);
                    };
                    let Some(&(r, i)) = by_map.get(&m) else {
                        return Some(4);
                    };
                    let r_atom = &reactants[r][i];
                    match (p_atom.tag, r_atom.tag) {
                        (Some(pt), Some(rt)) => {
                            let mut flag = if pt == rt { 2 } else { 1 };
                            if let Some(swaps) =
                                template_swaps(r_atom, &reactants[r], p_atom, template)
                                && swaps % 2 == 1
                            {
                                flag = 3 - flag;
                            }
                            Some(flag)
                        }
                        (Some(_), None) => Some(4),
                        (None, Some(_)) => Some(3),
                        (None, None) => None,
                    }
                })
                .collect()
        })
        .collect()
}

/// `countSwapsBetweenReactantAndProduct` on the two template atoms.
fn template_swaps(
    r_atom: &TemplateAtom,
    r_template: &[TemplateAtom],
    p_atom: &TemplateAtom,
    p_template: &[TemplateAtom],
) -> Option<usize> {
    let (rd, pd) = (r_atom.order.len(), p_atom.order.len());
    if rd < 3 || pd < 3 || rd.abs_diff(pd) > 1 {
        return None;
    }
    let neighbour_maps = |atom: &TemplateAtom, template: &[TemplateAtom], other: usize| {
        let mut order: Vec<i32> = atom
            .order
            .iter()
            .map(|&n| template[n].map.map_or(-1, i32::from))
            .collect();
        if atom.order.len() < other {
            order.push(-1);
        }
        let unmapped = order.iter().filter(|&&m| m < 0).count();
        (order, unmapped)
    };
    let (mut r_order, r_unmapped) = neighbour_maps(r_atom, r_template, pd);
    if r_unmapped > 1 {
        return None;
    }
    let (mut p_order, p_unmapped) = neighbour_maps(p_atom, p_template, rd);
    if p_unmapped > 1 {
        return None;
    }
    if !order_overlap(&mut r_order, r_unmapped, &p_order)
        || !order_overlap(&mut p_order, p_unmapped, &r_order)
    {
        return None;
    }
    count_swaps(&r_order, &p_order)
}

/// `checkOrderOverlap`.
fn order_overlap(order: &mut [i32], unmapped: usize, reference: &[i32]) -> bool {
    for &elem in reference {
        if elem >= 0 && !order.contains(&elem) {
            if unmapped == 0 {
                return false;
            }
            match order.iter().position(|&m| m == -1) {
                Some(slot) => order[slot] = elem,
                None => return false,
            }
        }
    }
    true
}

/// Where a product atom came from, for ordering its bonds.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum Node {
    /// A product-template atom.
    Template(usize),
    /// A reactant atom carried into the product.
    Carried(usize, AtomIdx),
}

/// RDKit's bond order of reactant atom `a`, as the SMILES parser created
/// it: the bond to the preceding atom, the bonds to atoms that follow it in
/// the text, then ring closures. Ring closures come last in RDKit but in
/// ring-number order, which a parsed molecule does not keep, so they are
/// returned separately for the caller to permute.
fn reactant_bond_order(mol: &Molecule, a: AtomIdx) -> (Vec<AtomIdx>, Vec<AtomIdx>) {
    let adjacency = mol.neighbor_slice(a);
    let mut prefix = Vec::with_capacity(adjacency.len());
    let mut rings = Vec::new();
    let parent = adjacency.first().map(|&(n, _)| n).filter(|&n| n < a);
    if let Some(p) = parent {
        prefix.push(p);
    }
    let mut children = Vec::new();
    for &(n, _) in adjacency {
        if Some(n) == parent {
            continue;
        }
        let is_child = n > a && mol.neighbor_slice(n).first().map(|&(p, _)| p) == Some(a);
        if is_child {
            children.push(n);
        } else {
            rings.push(n);
        }
    }
    children.sort_unstable();
    prefix.extend(children);
    (prefix, rings)
}

/// All orderings of `items` (Heap's algorithm), or `None` past the cap.
fn permutations<T: Clone>(items: &[T], cap: usize) -> Option<Vec<Vec<T>>> {
    let n = items.len();
    let count = (1..=n).try_fold(1usize, |acc, k| acc.checked_mul(k))?;
    if count > cap {
        return None;
    }
    let mut out = Vec::with_capacity(count);
    let mut a = items.to_vec();
    let mut c = vec![0usize; n];
    out.push(a.clone());
    let mut i = 0;
    while i < n {
        if c[i] < i {
            if i % 2 == 0 {
                a.swap(0, i);
            } else {
                a.swap(c[i], i);
            }
            out.push(a.clone());
            c[i] += 1;
            i = 0;
        } else {
            c[i] = 0;
            i += 1;
        }
    }
    Some(out)
}

/// The reactant atom's chiral tag referred to `order` (and an implicit H
/// last), from chematic's text-order stereo record.
fn reactant_tag(mol: &Molecule, a: AtomIdx, order: &[AtomIdx]) -> Option<bool> {
    let atom = mol.atom(a);
    let clockwise = match atom.chirality {
        Chirality::Clockwise => true,
        Chirality::CounterClockwise => false,
        _ => return None,
    };
    let recorded = mol.stereo_neighbor_order(a)?;
    let mut rd: Vec<u32> = order.iter().map(|n| n.0).collect();
    if recorded.contains(&STEREO_H_SENTINEL) {
        rd.push(STEREO_H_SENTINEL);
    }
    Some(clockwise == permutation_parity(recorded, &rd)?)
}

/// A chematic tetrahedral configuration: clockwise over a stereo neighbour
/// order (`None`: no tag).
type Config = Option<(bool, Vec<u32>)>;

/// One earlier traversal's bonds to a mapped atom: the traversal root's own
/// bond (if any), then carried atoms' bonds in an unknown order.
type BondGroup = (Option<Node>, Vec<Node>);

/// Bonds RDKit adds to one mapped atom while traversing the reactant from
/// mapped atoms that come earlier in the match: one group per earlier
/// atom, the atom's own bond (if any) first, then bonds from carried atoms
/// in an order that depends on the traversal.
#[derive(Default, Clone)]
struct Earlier {
    groups: Vec<BondGroup>,
    visited: Vec<bool>,
}

/// One reactant's atoms as the product sees them.
struct ReactantView {
    /// Reactant atom → product-template atom, for atoms mapped into it.
    mapped: FxHashMap<AtomIdx, usize>,
    /// Matched atoms not in this product.
    skipped: FxHashSet<AtomIdx>,
    /// Mapped atoms in match (template-atom) order.
    order: Vec<AtomIdx>,
}

/// Recompute the chiral tags of the product-template atoms of product `p`
/// as RDKit would, overriding the native tags on those atoms. Carried atoms
/// keep the native (geometric) stereo, which RDKit also preserves.
pub(crate) fn apply_product_stereo(
    info: &ReactionStereo,
    reactants: &[&Molecule],
    per_reactant: &[FxHashMap<usize, AtomIdx>],
    p: usize,
    template_idx_to_new: &[Option<AtomIdx>],
    src_to_new: &FxHashMap<(usize, AtomIdx), AtomIdx>,
    product: &mut Molecule,
) -> Result<(), StereoUnsupported> {
    let template = &info.products[p];
    let mut product_by_map: FxHashMap<u16, usize> = FxHashMap::default();
    let mut duplicate = false;
    for (i, atom) in template.iter().enumerate() {
        if let Some(m) = atom.map {
            duplicate |= product_by_map.insert(m, i).is_some();
        }
    }
    let mut source: Vec<Option<(usize, AtomIdx)>> = vec![None; template.len()];
    let mut views = Vec::with_capacity(reactants.len());
    for (k, matched) in per_reactant.iter().enumerate() {
        let mut pairs: Vec<(usize, AtomIdx)> = matched.iter().map(|(&t, &r)| (t, r)).collect();
        pairs.sort_unstable_by_key(|&(t, _)| t);
        let mut view = ReactantView {
            mapped: FxHashMap::default(),
            skipped: FxHashSet::default(),
            order: Vec::new(),
        };
        for (t, r) in pairs {
            let target = info.reactants[k]
                .get(t)
                .and_then(|a| a.map)
                .and_then(|m| product_by_map.get(&m).copied());
            match target {
                Some(i) => {
                    duplicate |= source[i].replace((k, r)).is_some();
                    view.mapped.insert(r, i);
                    view.order.push(r);
                }
                None => {
                    view.skipped.insert(r);
                }
            }
        }
        views.push(view);
    }

    // Carried atoms keep their geometry (`checkAndCorrectChiralityOfProduct`):
    // re-express the reactant's stereo record in product atoms so it does
    // not depend on the product's bond order.
    for (&(k, r), &idx) in src_to_new {
        if template_idx_to_new.contains(&Some(idx)) {
            continue;
        }
        let atom = reactants[k].atom(r);
        if !atom.chirality.is_tetrahedral() {
            continue;
        }
        let remapped: Option<Vec<u32>> = reactants[k].stereo_neighbor_order(r).and_then(|order| {
            order
                .iter()
                .map(|&n| {
                    if n == STEREO_H_SENTINEL {
                        Some(n)
                    } else {
                        src_to_new.get(&(k, AtomIdx(n))).map(|a| a.0)
                    }
                })
                .collect()
        });
        let degree_kept = reactants[k].degree(r) == product.degree(idx);
        match remapped {
            Some(order) if degree_kept => {
                product.set_chirality(idx, atom.chirality);
                product.set_stereo_neighbor_order(idx, order);
            }
            _ => product.set_chirality(idx, Chirality::None),
        }
    }

    let chiral_source =
        |i: usize| source[i].is_some_and(|(k, r)| reactants[k].atom(r).chirality.is_tetrahedral());
    let involved: Vec<usize> = (0..template.len())
        .filter(|&i| chiral_source(i) || template[i].tag.is_some())
        .collect();
    if involved.is_empty() {
        // Neither a template tag nor a chiral mapped atom: RDKit leaves
        // every template atom without a tag, as does the native build.
        for &idx in template_idx_to_new.iter().flatten() {
            if product.atom(idx).chirality != Chirality::None {
                product.set_chirality(idx, Chirality::None);
            }
        }
        return Ok(());
    }
    if duplicate {
        return Err(StereoUnsupported::BondOrder);
    }

    let template_bond = |i: usize, j: usize| template[i].order.contains(&j);
    let mut earlier: FxHashMap<(usize, AtomIdx), Earlier> = FxHashMap::default();
    for (k, view) in views.iter().enumerate() {
        if !view.order.iter().any(|r| {
            let i = view.mapped[r];
            involved.contains(&i)
        }) {
            continue;
        }
        let mol = reactants[k];
        let mut visited = vec![false; mol.atom_count()];
        let mut added = vec![false; mol.atom_count()];
        let mut pending: FxHashMap<AtomIdx, Vec<BondGroup>> = FxHashMap::default();
        for &m in &view.order {
            let mi = view.mapped[&m];
            earlier.insert(
                (k, m),
                Earlier {
                    groups: pending.remove(&m).unwrap_or_default(),
                    visited: visited.clone(),
                },
            );
            // `addReactantNeighborsToProduct` from `m`: the bonds it adds to
            // mapped atoms whose own traversal comes later.
            let mut round: FxHashMap<AtomIdx, (Option<Node>, Vec<Node>)> = FxHashMap::default();
            let mut queue = std::collections::VecDeque::from([m]);
            while let Some(l) = queue.pop_front() {
                visited[l.0 as usize] = true;
                for (n, _) in mol.neighbors(l) {
                    let ni = n.0 as usize;
                    if visited[ni] || view.skipped.contains(&n) {
                        continue;
                    }
                    if let Some(&nj) = view.mapped.get(&n) {
                        let group = round.entry(n).or_default();
                        if l == m {
                            // Both mapped (github #1387): bonded unless the
                            // templates already decide this pair.
                            let ma = template[mi].map.unwrap_or(0);
                            let na = template[nj].map.unwrap_or(0);
                            if !template_bond(mi, nj)
                                && !info.reactant_bond_maps.contains(&(ma.min(na), ma.max(na)))
                            {
                                group.0 = Some(Node::Template(mi));
                            }
                        } else {
                            group.1.push(Node::Carried(k, l));
                        }
                    } else if !added[ni] {
                        added[ni] = true;
                        queue.push_back(n);
                    }
                }
            }
            for (n, group) in round {
                if group.0.is_some() || !group.1.is_empty() {
                    pending.entry(n).or_default().push(group);
                }
            }
        }
    }

    // Per involved atom: the RDKit result under every bond-order option.
    let mut decided: Vec<(AtomIdx, Vec<Config>)> = Vec::new();
    for i in 0..template.len() {
        let Some(new_idx) = template_idx_to_new.get(i).copied().flatten() else {
            continue;
        };
        let flag = info.flags[p][i];
        let template_part: Vec<Node> = template[i]
            .order
            .iter()
            .map(|&j| Node::Template(j))
            .collect();
        let results: Vec<Option<(bool, Vec<Node>)>> = match source[i] {
            None => vec![
                (flag == Some(4))
                    .then_some(template[i].tag)
                    .flatten()
                    .map(|t| (t, template_part.clone())),
            ],
            Some((k, a)) => {
                let mol = reactants[k];
                let chiral = mol.atom(a).chirality.is_tetrahedral();
                let created = flag == Some(4) && template[i].tag.is_some();
                if !chiral && !created {
                    vec![None]
                } else {
                    let view = &views[k];
                    let e = earlier.get(&(k, a)).cloned().unwrap_or_default();
                    let (prefix, rings) = reactant_bond_order(mol, a);
                    let ring_orders = permutations(&rings, MAX_ORDER_OPTIONS)
                        .ok_or(StereoUnsupported::BondOrder)?;
                    let mut group_options: Vec<Vec<Node>> = vec![Vec::new()];
                    for (fixed, carried) in &e.groups {
                        let perms = permutations(carried, MAX_ORDER_OPTIONS)
                            .ok_or(StereoUnsupported::BondOrder)?;
                        let mut next = Vec::new();
                        for base in &group_options {
                            for perm in &perms {
                                let mut v = base.clone();
                                v.extend(fixed.iter().copied());
                                v.extend(perm.iter().copied());
                                next.push(v);
                            }
                        }
                        if next.len() * ring_orders.len() > MAX_ORDER_OPTIONS {
                            return Err(StereoUnsupported::BondOrder);
                        }
                        group_options = next;
                    }
                    let mut out = Vec::new();
                    for ring_order in &ring_orders {
                        let mut r_order = prefix.clone();
                        r_order.extend(ring_order.iter().copied());
                        let own: Vec<Node> = r_order
                            .iter()
                            .filter(|n| {
                                !e.visited.get(n.0 as usize).copied().unwrap_or(false)
                                    && !view.skipped.contains(n)
                            })
                            .filter_map(|&n| match view.mapped.get(&n) {
                                Some(&j) => {
                                    let ma = template[i].map.unwrap_or(0);
                                    let na = template[j].map.unwrap_or(0);
                                    (!template_bond(i, j)
                                        && !info
                                            .reactant_bond_maps
                                            .contains(&(ma.min(na), ma.max(na))))
                                    .then_some(Node::Template(j))
                                }
                                None => Some(Node::Carried(k, n)),
                            })
                            .collect();
                        for groups in &group_options {
                            let mut list = template_part.clone();
                            list.extend(groups.iter().copied());
                            list.extend(own.iter().copied());
                            out.push(rdkit_tag(
                                info, p, i, flag, mol, k, a, &r_order, &list, &source,
                            )?);
                        }
                    }
                    out
                }
            }
        };
        // The distinct answers over all options (normally one).
        let mut outcomes: Vec<(Option<bool>, Config)> = Vec::new();
        for result in results {
            let converted = match result {
                None => None,
                Some((cw, list)) => Some(to_product_config(
                    cw,
                    &list,
                    new_idx,
                    product,
                    template_idx_to_new,
                    src_to_new,
                )?),
            };
            let key = converted.as_ref().map(|(cw, order)| {
                let mut sorted = order.clone();
                sorted.sort_unstable();
                *cw == permutation_parity(order, &sorted).unwrap_or(true)
            });
            if !outcomes.iter().any(|(k, _)| *k == key) {
                outcomes.push((key, converted));
            }
        }
        decided.push((new_idx, outcomes.into_iter().map(|(_, c)| c).collect()));
    }
    let ambiguous: Vec<usize> = (0..decided.len())
        .filter(|&d| decided[d].1.len() > 1)
        .collect();
    for (idx, outcomes) in &decided {
        set_config(product, *idx, &outcomes[0]);
    }
    if ambiguous.is_empty() {
        return Ok(());
    }
    // Some atoms' tags depend on an unrecorded bond order. Accept when the
    // product is the same molecule whichever way they fall (the atom is not
    // a stereocentre in the product), else decline.
    let combinations = ambiguous
        .iter()
        .try_fold(1usize, |n, &d| n.checked_mul(decided[d].1.len()))
        .filter(|&n| n <= 64)
        .ok_or(StereoUnsupported::BondOrder)?;
    let reference = chematic_smiles::canonical_smiles(product);
    for combination in 1..combinations {
        let mut variant = product.clone();
        let mut rest = combination;
        for &d in &ambiguous {
            let (idx, outcomes) = &decided[d];
            set_config(&mut variant, *idx, &outcomes[rest % outcomes.len()]);
            rest /= outcomes.len();
        }
        if chematic_smiles::canonical_smiles(&variant) != reference {
            return Err(StereoUnsupported::BondOrder);
        }
    }
    Ok(())
}

fn set_config(product: &mut Molecule, idx: AtomIdx, config: &Config) {
    match config {
        Some((cw, order)) => {
            product.set_chirality(
                idx,
                if *cw {
                    Chirality::Clockwise
                } else {
                    Chirality::CounterClockwise
                },
            );
            product.set_stereo_neighbor_order(idx, order.clone());
        }
        None => {
            if product.atom(idx).chirality != Chirality::None {
                product.set_chirality(idx, Chirality::None);
            }
        }
    }
}

/// The chiral tag RDKit leaves on product atom `i` with bond list `list`.
#[allow(clippy::too_many_arguments)]
fn rdkit_tag(
    info: &ReactionStereo,
    p: usize,
    i: usize,
    flag: Option<u8>,
    mol: &Molecule,
    k: usize,
    a: AtomIdx,
    r_order: &[AtomIdx],
    list: &[Node],
    source: &[Option<(usize, AtomIdx)>],
) -> Result<Option<(bool, Vec<Node>)>, StereoUnsupported> {
    let template = &info.products[p];
    // `convertTemplateToMol`: only flag 4 keeps the template's tag.
    let mut tag = if flag == Some(4) {
        template[i].tag
    } else {
        None
    };
    let reactant = if mol.atom(a).chirality.is_tetrahedral() {
        Some(reactant_tag(mol, a, r_order).ok_or(StereoUnsupported::BondOrder)?)
    } else {
        None
    };
    let Some(rt) = reactant else {
        return Ok(tag.map(|t| (t, list.to_vec())));
    };
    // `checkProductChirality`.
    match flag {
        Some(1) => tag = Some(!rt),
        Some(2) => tag = Some(rt),
        Some(3) => tag = None,
        _ => {}
    }
    // `checkAndCorrectChiralityOfMatchingAtomsInProduct`.
    if !matches!(flag, Some(3 | 4)) {
        let (rdeg, pdeg) = (r_order.len(), list.len());
        if rdeg >= 3 && pdeg >= 3 && rdeg.abs_diff(pdeg) <= 1 {
            let mut unknown = 0;
            let mut p_order: Vec<Option<AtomIdx>> = Vec::with_capacity(pdeg);
            for node in list {
                let traced = match *node {
                    Node::Template(j) => source[j].filter(|&(kk, _)| kk == k).map(|(_, b)| b),
                    Node::Carried(kk, b) => (kk == k).then_some(b),
                }
                .filter(|&b| r_order.contains(&b));
                if traced.is_none() {
                    unknown += 1;
                    if unknown > 1 {
                        break;
                    }
                }
                p_order.push(traced);
            }
            if unknown == 1 {
                if rdeg == pdeg {
                    let missing = r_order
                        .iter()
                        .copied()
                        .find(|b| !p_order.contains(&Some(*b)));
                    if let (Some(m), Some(slot)) =
                        (missing, p_order.iter().position(Option::is_none))
                    {
                        p_order[slot] = Some(m);
                    }
                    unknown = 0;
                } else if pdeg > rdeg {
                    p_order.retain(Option::is_some);
                    unknown = 0;
                }
            }
            if unknown == 0 {
                let p_order: Vec<AtomIdx> = p_order.into_iter().flatten().collect();
                let mut r = r_order.to_vec();
                if rdeg > pdeg {
                    r.retain(|b| p_order.contains(b));
                }
                let swaps = count_swaps(&r, &p_order).ok_or(StereoUnsupported::BondOrder)?;
                let mut t = rt;
                if swaps % 2 == 1 {
                    t = !t;
                }
                if flag == Some(1) {
                    t = !t;
                }
                tag = Some(t);
            }
        }
    }
    Ok(tag.map(|t| (t, list.to_vec())))
}

/// Express an RDKit tag over `list` (an implicit H last) as chematic's
/// chirality and stereo neighbour order on product atom `idx`.
fn to_product_config(
    clockwise: bool,
    list: &[Node],
    idx: AtomIdx,
    product: &Molecule,
    template_idx_to_new: &[Option<AtomIdx>],
    src_to_new: &FxHashMap<(usize, AtomIdx), AtomIdx>,
) -> Result<(bool, Vec<u32>), StereoUnsupported> {
    let mut order = Vec::with_capacity(list.len() + 1);
    for node in list {
        let n = match *node {
            Node::Template(j) => template_idx_to_new.get(j).copied().flatten(),
            Node::Carried(k, b) => src_to_new.get(&(k, b)).copied(),
        }
        .ok_or(StereoUnsupported::BondOrder)?;
        order.push(n.0);
    }
    let actual: FxHashSet<u32> = product.neighbors(idx).map(|(n, _)| n.0).collect();
    let listed: FxHashSet<u32> = order.iter().copied().collect();
    if actual != listed || listed.len() != order.len() {
        return Err(StereoUnsupported::BondOrder);
    }
    if order.len() == 3 && chematic_core::implicit_hcount(product, idx) > 0 {
        order.push(STEREO_H_SENTINEL);
    }
    Ok((clockwise, order))
}

/// The reactant as RDKit's SMILES parser leaves it: tetrahedral tags on
/// atoms that cannot be stereocentres are removed (legacy stereo perception
/// with `cleanIt`, the 2026.03.6 default), together with a bracket H written
/// only for the tag. `None` when nothing changes.
///
/// The degree and hydrogen rules follow `isAtomPotentialChiralCenter`. Where
/// RDKit compares CIP ranks refined by already assigned centres, this asks
/// whether inverting the tag gives a different molecule (by stereo-aware
/// canonical SMILES), which covers ring cis/trans pairs and
/// pseudo-asymmetric centres the same way.
pub(crate) fn rdkit_parse_cleanup(mol: &Molecule) -> Option<Molecule> {
    let tagged: Vec<AtomIdx> = mol
        .atoms()
        .filter(|(_, a)| a.chirality.is_tetrahedral())
        .map(|(i, _)| i)
        .collect();
    if tagged.is_empty() {
        return None;
    }
    let classes = chematic_smiles::topological_equivalence_classes(mol);
    let legal = |a: AtomIdx| -> bool {
        let atom = mol.atom(a);
        let z = atom.element.atomic_number();
        let degree = mol.degree(a);
        let h = chematic_core::implicit_hcount(mol, a) as usize;
        if degree + h > 4 || degree + h < 3 {
            false
        } else if degree < 3 {
            matches!(z, 15 | 33)
        } else if degree == 3 && h != 1 {
            match z {
                7 => {
                    let nbrs: Vec<AtomIdx> = mol.neighbors(a).map(|(n, _)| n).collect();
                    let three_ring = nbrs.iter().enumerate().any(|(i, &x)| {
                        nbrs[i + 1..]
                            .iter()
                            .any(|&y| mol.neighbors(x).any(|(n, _)| n == y))
                    });
                    three_ring
                        && mol
                            .neighbors(a)
                            .all(|(_, b)| mol.bond(b).order == chematic_core::BondOrder::Single)
                }
                15 | 33 => true,
                16 | 34 => {
                    let valence: u32 = mol
                        .neighbors(a)
                        .map(|(_, b)| u32::from(mol.bond(b).order.order_int()))
                        .sum();
                    valence == 4 || (valence == 3 && atom.charge == 1)
                }
                _ => false,
            }
        } else {
            true
        }
    };
    let duplicate_neighbours = |a: AtomIdx| {
        let mut seen = FxHashSet::default();
        !mol.neighbors(a)
            .all(|(n, _)| seen.insert(classes[n.0 as usize]))
    };
    let mut out = mol.clone();
    let mut changed = false;
    let clear = |out: &mut Molecule, a: AtomIdx| {
        out.set_chirality(a, Chirality::None);
        // Issue 194: the bracket H was there only for the tag.
        let atom = mol.atom(a);
        if atom.hydrogen_count == Some(1) && atom.charge == 0 && !atom.aromatic {
            out.set_hydrogen_count(a, None);
        }
    };
    let mut undecided = Vec::new();
    for &a in &tagged {
        if !legal(a) {
            clear(&mut out, a);
            changed = true;
        } else if duplicate_neighbours(a) {
            undecided.push(a);
        }
    }
    // RDKit's ring special case (`findChiralAtomSpecialCases`): tagged
    // candidates sharing a ring system keep their tags together, as in
    // 1,4-disubstituted cyclohexanes, even where one of them is redundant.
    let ring_bond = ring_bonds(mol);
    let candidate = |a: AtomIdx| -> bool {
        let atom = mol.atom(a);
        if atom.element.atomic_number() == 7
            && mol.degree(a) + chematic_core::implicit_hcount(mol, a) as usize == 3
        {
            return false;
        }
        let mut ring_nbrs = Vec::new();
        let mut other = Vec::new();
        for (n, b) in mol.neighbors(a) {
            if ring_bond[b.0 as usize] {
                ring_nbrs.push(classes[n.0 as usize]);
            } else {
                other.push(classes[n.0 as usize]);
            }
        }
        let mut distinct = ring_nbrs.clone();
        distinct.sort_unstable();
        distinct.dedup();
        match other.len() {
            _ if ring_nbrs.is_empty() => false,
            2 => other[0] != other[1] && ring_nbrs.len() != distinct.len(),
            1 => ring_nbrs.len() > distinct.len(),
            0 => {
                (ring_nbrs.len() == 4 && distinct.len() == 3)
                    || (ring_nbrs.len() == 3 && distinct.len() == 2)
            }
            _ => false,
        }
    };
    let special: Vec<AtomIdx> = undecided
        .iter()
        .copied()
        .filter(|&a| candidate(a))
        .collect();
    let mut paired = FxHashSet::default();
    for &a in &special {
        let mut seen = FxHashSet::from_iter([a]);
        let mut queue = vec![a];
        while let Some(x) = queue.pop() {
            for (n, b) in mol.neighbors(x) {
                if ring_bond[b.0 as usize] && seen.insert(n) {
                    queue.push(n);
                }
            }
        }
        if special.iter().any(|&o| o != a && seen.contains(&o)) {
            paired.insert(a);
        }
    }
    undecided.retain(|a| !paired.contains(a));
    // A tag on an atom with two alike neighbours matters only if inverting
    // it gives another molecule; repeat while clearing changes the answer.
    loop {
        let mut cleared_any = false;
        let reference = chematic_smiles::canonical_smiles(&out);
        undecided.retain(|&a| {
            let mut inverted = out.clone();
            let flipped = match out.atom(a).chirality {
                Chirality::Clockwise => Chirality::CounterClockwise,
                _ => Chirality::Clockwise,
            };
            inverted.set_chirality(a, flipped);
            if chematic_smiles::canonical_smiles(&inverted) == reference {
                clear(&mut out, a);
                cleared_any = true;
                false
            } else {
                true
            }
        });
        if !cleared_any {
            break;
        }
        changed = true;
    }
    changed.then_some(out)
}

/// Whether each bond lies on a ring (is not a bridge).
fn ring_bonds(mol: &Molecule) -> Vec<bool> {
    let n = mol.atom_count();
    let mut order = vec![usize::MAX; n];
    let mut low = vec![0usize; n];
    let mut ring = vec![true; mol.bond_count()];
    let mut counter = 0;
    for root in 0..n {
        if order[root] != usize::MAX {
            continue;
        }
        let mut stack: Vec<(usize, Option<usize>, usize)> = vec![(root, None, 0)];
        order[root] = counter;
        low[root] = counter;
        counter += 1;
        while let Some(top) = stack.last_mut() {
            let (v, via, k) = *top;
            if let Some(&(w, b)) = mol.neighbor_slice(AtomIdx(v as u32)).get(k) {
                top.2 += 1;
                let (w, b) = (w.0 as usize, b.0 as usize);
                if Some(b) == via {
                    continue;
                }
                if order[w] == usize::MAX {
                    order[w] = counter;
                    low[w] = counter;
                    counter += 1;
                    stack.push((w, Some(b), 0));
                } else {
                    low[v] = low[v].min(order[w]);
                }
            } else {
                stack.pop();
                if let Some(&(p, _, _)) = stack.last() {
                    low[p] = low[p].min(low[v]);
                    if low[v] > order[p] {
                        ring[via.expect("tree edge")] = false;
                    }
                }
            }
        }
    }
    ring
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tokenizer_orders_ring_closures_last_by_number() {
        let atoms = tokenize_template("C[C@@H]1CC[C@H](O)CC1").unwrap();
        assert_eq!(atoms[1].order, vec![0, 2, 7]);
        assert_eq!(atoms[1].tag, Some(false));
        assert_eq!(atoms[4].order, vec![3, 5, 6]);
        assert_eq!(atoms[4].tag, Some(false));
    }

    #[test]
    fn tokenizer_h_only_in_simple_spelling() {
        assert_eq!(
            tokenize_template("[C@@H](N)(C)O").unwrap()[0].tag,
            Some(false)
        );
        assert_eq!(
            tokenize_template("[C@@;H1:2](N)(C)O").unwrap()[0].tag,
            Some(true)
        );
        assert_eq!(
            tokenize_template("N[C@@H:2](C)O").unwrap()[1].tag,
            Some(true)
        );
        assert_eq!(
            tokenize_template("[C@@H:2](N)(C)O").unwrap()[0].map,
            Some(2)
        );
    }

    #[test]
    fn count_swaps_matches_rdkit() {
        assert_eq!(count_swaps(&[0, 1, 2, 3], &[1, 0, 2, 3]), Some(1));
        assert_eq!(count_swaps(&[0, 1, 2, 3], &[1, 2, 3, 0]), Some(3));
        assert_eq!(count_swaps(&[0, 1, 2, 3], &[1, 2, 0, 3]), Some(2));
    }
}
