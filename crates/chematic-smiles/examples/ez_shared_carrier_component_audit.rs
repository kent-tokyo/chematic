//! Wave 2B (issue #149) pre-implementation audit: measures the shared-
//! candidate-bond coupling graph BEFORE any solver is designed, so the
//! design is grounded in the actual current data (this repo's molecules),
//! not just the issue text's own description.
//!
//! A "stereo alkene end" node is a double-bond terminus with exactly 2
//! non-double-bond substituents -- the same ambiguity precondition
//! `chematic_smiles::canonical`'s private `resolve_ez_marker_for_end` uses
//! (`subs.len() == 2`). Re-derived here using only `chematic_core`'s public
//! `Molecule` API (bonds/neighbors/bond order) -- a purely topological
//! survey needs no access to `canonical.rs`'s private `CanonicalWriter`.
//!
//! An edge connects two such nodes whenever they are directly bonded via a
//! non-double bond: that physical bond is then a candidate substituent
//! (carrier) for BOTH ends' resolution at once (the only way a single bond
//! can be shared between two independently-stereogenic double bonds' own
//! marker-carrier choices). Since each node has at most 2 candidate
//! substituent bonds, every node has degree <= 2 in this graph -- so every
//! connected component is a simple path or a simple cycle. Reports the
//! size distribution (of components with >=2 nodes -- singletons have no
//! coupling to resolve) plus, per component, whether it is a path/cycle and
//! whether any interior node has zero "private" (non-shared) substituents.
//!
//! Run against the 18 pinned residual fixtures (no argument) or a full
//! corpus (any one-SMILES-per-line file -- `scripts/descriptor_census_
//! corpus.smi`, committed to this repo, needs no external download):
//!
//! ```text
//! cargo run -p chematic-smiles --release --example ez_shared_carrier_component_audit
//! cargo run -p chematic-smiles --release --example ez_shared_carrier_component_audit -- scripts/descriptor_census_corpus.smi
//! ```

use std::collections::{HashMap, HashSet};
use std::env;
use std::fs;

use chematic_core::{AtomIdx, BondIdx, BondOrder, Molecule};

const EZ_SHARED_CANDIDATE_BOND_RESIDUALS: &[&str] = &[
    r"CCCCC/N=c1\c(O)c(O)\c1=N/[C@@H](Cc1ccc(NC(=O)c2c(Cl)cncc2Cl)cc1)C(=O)O",
    r"O=C(Nc1ccc(C[C@H](/N=c2\c(O)c(O)\c2=N/Cc2ccccc2)C(=O)O)cc1)c1c(Cl)cncc1Cl",
    r"CCC/N=c1\c(O)c(O)\c1=N/[C@@H](Cc1ccc(NC(=O)c2c(Cl)cncc2Cl)cc1)C(=O)O",
    r"O=C(Nc1ccc(C[C@H](/N=c2\c(O)c(O)\c2=N/c2ccccc2)C(=O)O)cc1)c1c(Cl)cncc1Cl",
    r"CC(C)(C)/N=c1\c(O)c(O)\c1=N/[C@@H](Cc1ccc(NC(=O)c2c(Cl)cncc2Cl)cc1)C(=O)O",
    r"CC1=C2CC[C@H](/C=N/N=C(N)N)[C@@]2(C)CC/C1=N\N=C(N)N",
    r"CC1=C2CC[C@@H](/C=N/N=C(N)N)[C@@]2(C)CC/C1=N\N=C(N)N",
    r"COC(=O)/C=C/[C@H]1CCC2=C(C)/C(=N/N=C(N)N)CC[C@@]21C",
    r"CCOC(=O)C1=C(C)N=C(C)/C(=C(/O)OCC)C1c1cccc(I)c1",
    r"CCOC(=O)C1=C(C)N=C(C)/C(=C(/O)OCC)C1c1ccc(I)cc1",
    r"CCO/C(O)=C1\C(C)=NC(C)=C(C(=O)OC)C1c1ccccc1C(F)(F)F",
    r"CCOC(=O)C1=C(C)N=C(C)/C(=C(/O)OCC)C1c1ccccc1OC",
    r"CCO/C(O)=C1\C(C)=NC(C)=C(C(=O)OC)C1c1ccccc1[N+](=O)[O-]",
    r"CCO/C(O)=C1\C(C)=NC(C)=C(C(=O)OC)C1c1ccc([N+](=O)[O-])cc1",
    r"CCOC(=O)C1=C(C)N=C(C)/C(=C(/O)OCC)C1c1cccc(C)c1",
    r"CCOC(=O)C1=C(C)N=C(C)/C(=C(/O)OCC)C1c1cccc(OC)c1",
    r"CCO/C(O)=C1\C(C)=NC(C)=C(C(=O)OC)C1c1cccc(C(F)(F)F)c1",
    r"CCO/C(O)=C(\C1=NCCN1)c1nnc(N)s1",
];

fn end_has_substituent(mol: &Molecule, end: AtomIdx) -> bool {
    mol.neighbors(end)
        .any(|(_, b)| mol.bond(b).order != BondOrder::Double)
}

fn substituents(mol: &Molecule, end: AtomIdx) -> Vec<(AtomIdx, BondIdx)> {
    mol.neighbors(end)
        .filter(|&(_, b)| mol.bond(b).order != BondOrder::Double)
        .collect()
}

/// Every atom that is a stereogenic double bond's terminus with exactly 2
/// candidate substituents -- mirrors `canonical.rs`'s private
/// `stereo_alkene_ends` computation exactly (both ends of the double bond
/// must have >=1 substituent for either end to count at all).
fn stereo_alkene_end_nodes(mol: &Molecule) -> HashSet<AtomIdx> {
    let mut nodes = HashSet::new();
    for (_, bond) in mol.bonds() {
        if bond.order != BondOrder::Double {
            continue;
        }
        if !end_has_substituent(mol, bond.atom1) || !end_has_substituent(mol, bond.atom2) {
            continue;
        }
        for end in [bond.atom1, bond.atom2] {
            if substituents(mol, end).len() == 2 {
                nodes.insert(end);
            }
        }
    }
    nodes
}

struct ComponentReport {
    size: usize,
    is_cycle: bool,
    any_node_zero_private: bool,
}

/// Connected components of the shared-candidate-bond coupling graph
/// (nodes = `stereo_alkene_end_nodes`, edges = direct non-double bonds
/// between two nodes). Only components with >=2 nodes are coupling at all.
fn coupling_components(mol: &Molecule) -> Vec<ComponentReport> {
    let nodes = stereo_alkene_end_nodes(mol);
    if nodes.len() < 2 {
        return Vec::new();
    }

    let mut adjacency: HashMap<AtomIdx, Vec<AtomIdx>> = HashMap::new();
    for &n in &nodes {
        let subs = substituents(mol, n);
        for &(sub_atom, _) in &subs {
            if nodes.contains(&sub_atom) {
                adjacency.entry(n).or_default().push(sub_atom);
            }
        }
    }

    let mut visited: HashSet<AtomIdx> = HashSet::new();
    let mut reports = Vec::new();
    for &start in &nodes {
        if visited.contains(&start) {
            continue;
        }
        // BFS over the coupling graph.
        let mut queue = vec![start];
        let mut component = Vec::new();
        visited.insert(start);
        while let Some(cur) = queue.pop() {
            component.push(cur);
            for &nb in adjacency.get(&cur).map(|v| v.as_slice()).unwrap_or(&[]) {
                if visited.insert(nb) {
                    queue.push(nb);
                }
            }
        }
        if component.len() < 2 {
            continue; // isolated node with no coupling edge at all
        }
        let n_edges: usize = component
            .iter()
            .map(|c| adjacency.get(c).map(|v| v.len()).unwrap_or(0))
            .sum::<usize>()
            / 2;
        let is_cycle = n_edges == component.len();
        let any_node_zero_private = component.iter().any(|&n| {
            let subs = substituents(mol, n);
            subs.iter().all(|&(a, _)| nodes.contains(&a))
        });
        reports.push(ComponentReport {
            size: component.len(),
            is_cycle,
            any_node_zero_private,
        });
    }
    reports
}

fn analyze(label: &str, smiles_list: &[String]) {
    let mut size_histogram: HashMap<usize, usize> = HashMap::new();
    let mut n_cycles = 0usize;
    let mut n_zero_private_components = 0usize;
    let mut n_molecules_with_coupling = 0usize;
    let mut parse_failures = 0usize;
    let mut max_component_size = 0usize;

    for s in smiles_list {
        let mol = match chematic_smiles::parse(s) {
            Ok(m) => m,
            Err(_) => {
                parse_failures += 1;
                continue;
            }
        };
        let components = coupling_components(&mol);
        if !components.is_empty() {
            n_molecules_with_coupling += 1;
        }
        for c in &components {
            *size_histogram.entry(c.size).or_insert(0) += 1;
            max_component_size = max_component_size.max(c.size);
            if c.is_cycle {
                n_cycles += 1;
            }
            if c.any_node_zero_private {
                n_zero_private_components += 1;
            }
        }
    }

    println!("=== {label} ===");
    println!("  molecules scanned: {}", smiles_list.len());
    println!("  parse failures: {parse_failures}");
    println!("  molecules with >=1 coupling component (size>=2): {n_molecules_with_coupling}");
    println!("  max component size observed: {max_component_size}");
    let mut sizes: Vec<_> = size_histogram.into_iter().collect();
    sizes.sort_by_key(|&(k, _)| k);
    println!("  component-size histogram (size -> count): {sizes:?}");
    println!("  components that are cycles (not paths): {n_cycles}");
    println!(
        "  components with >=1 interior node having ZERO private substituents: {n_zero_private_components}"
    );
}

fn main() {
    let fixtures: Vec<String> = EZ_SHARED_CANDIDATE_BOND_RESIDUALS
        .iter()
        .map(|s| s.to_string())
        .collect();
    analyze(
        "18 pinned EZ_SHARED_CANDIDATE_BOND_RESIDUALS fixtures",
        &fixtures,
    );

    if let Some(path) = env::args().nth(1) {
        let text = fs::read_to_string(&path).unwrap_or_else(|e| panic!("read {path}: {e}"));
        let corpus: Vec<String> = text
            .lines()
            .map(|l| l.split_whitespace().next().unwrap_or("").to_string())
            .filter(|l| !l.is_empty())
            .collect();
        analyze(&format!("full corpus ({path})"), &corpus);
    } else {
        eprintln!(
            "(no corpus path given -- pass one as an argument for the full-corpus histogram)"
        );
    }
}
