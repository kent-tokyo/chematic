//! Exports `kekulize` inputs and results for the Lean differential checker.
//!
//! ```text
//! kekule-diff-harness random <seed> <count>   # random graphs, arbitrary must-sets
//! kekule-diff-harness smiles <file>...        # real molecules, real chemistry rules
//! ```
//!
//! Output format is documented in `../lean/Main.lean`.

use std::io::{BufWriter, Write};

use chematic_core::{
    Atom, BondOrder, Element, Molecule, atom_must_be_matched, kekulize, kekulize_with,
    molecule::{AtomIdx, MoleculeBuilder},
};

/// Small deterministic PRNG (SplitMix64) so runs are reproducible without deps.
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
    fn chance(&mut self, num: u64, den: u64) -> bool {
        self.below(den) < num
    }
}

/// Random edge set over `n` atoms: either G(n, p) or a fused-ring system
/// (rings sharing an edge, plus occasional chords), which is where odd cycles
/// and blossoms show up.
fn random_edges(rng: &mut Rng) -> (usize, Vec<(usize, usize)>) {
    let mut edges: Vec<(usize, usize)> = Vec::new();
    let add = |edges: &mut Vec<(usize, usize)>, a: usize, b: usize| {
        if a != b && !edges.iter().any(|&(x, y)| (x, y) == (a, b) || (x, y) == (b, a)) {
            edges.push((a, b));
        }
    };
    if rng.chance(1, 3) {
        let n = 2 + rng.below(13) as usize;
        let p = 15 + rng.below(40);
        for a in 0..n {
            for b in a + 1..n {
                if rng.below(100) < p {
                    add(&mut edges, a, b);
                }
            }
        }
        return (n, edges);
    }
    // Fused rings.
    let first = 3 + rng.below(6) as usize;
    let mut n = first;
    for i in 0..first {
        add(&mut edges, i, (i + 1) % first);
    }
    for _ in 0..rng.below(4) {
        // Fuse a new ring of size k onto a random existing edge (a, b).
        let (a, b) = edges[rng.below(edges.len() as u64) as usize];
        let k = 3 + rng.below(6) as usize;
        let mut prev = a;
        for _ in 0..k - 2 {
            add(&mut edges, prev, n);
            prev = n;
            n += 1;
        }
        add(&mut edges, prev, b);
    }
    for _ in 0..rng.below(3) {
        if rng.chance(1, 3) {
            let a = rng.below(n as u64) as usize;
            let b = rng.below(n as u64) as usize;
            add(&mut edges, a, b);
        }
    }
    (n, edges)
}

fn emit_case(
    out: &mut impl Write,
    id: &str,
    mol: &Molecule,
    must: &[bool],
    excl: &[u32],
    result: Result<chematic_core::KekuleResult, chematic_core::KekuleError>,
) -> std::io::Result<()> {
    let aromatic: Vec<_> = mol
        .bonds()
        .filter(|(_, b)| b.order == BondOrder::Aromatic)
        .collect();
    writeln!(out, "c {id}")?;
    write!(out, "e")?;
    for (_, b) in &aromatic {
        write!(out, " {} {}", b.atom1.0, b.atom2.0)?;
    }
    writeln!(out)?;
    write!(out, "m")?;
    for (i, _) in must.iter().enumerate().filter(|(_, m)| **m) {
        write!(out, " {i}")?;
    }
    writeln!(out)?;
    write!(out, "x")?;
    for x in excl {
        write!(out, " {x}")?;
    }
    writeln!(out)?;
    match result {
        Ok(map) => {
            write!(out, "r ok")?;
            for (k, (bidx, _)) in aromatic.iter().enumerate() {
                if map.get(bidx) == Some(&BondOrder::Double) {
                    write!(out, " {k}")?;
                }
            }
            writeln!(out)?;
        }
        Err(_) => writeln!(out, "r err")?,
    }
    writeln!(out, "end")
}

fn aromatic_atoms(mol: &Molecule) -> Vec<bool> {
    let mut on = vec![false; mol.atom_count()];
    for (_, b) in mol.bonds() {
        if b.order == BondOrder::Aromatic {
            on[b.atom1.0 as usize] = true;
            on[b.atom2.0 as usize] = true;
        }
    }
    on
}

fn random_mode(seed: u64, count: usize, out: &mut impl Write) -> std::io::Result<()> {
    let mut rng = Rng(seed);
    for case in 0..count {
        let (n, edges) = random_edges(&mut rng);
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<AtomIdx> = (0..n).map(|_| b.add_atom(Atom::aromatic(Element::C))).collect();
        for &(x, y) in &edges {
            b.add_bond(atoms[x], atoms[y], BondOrder::Aromatic).unwrap();
        }
        let mol = b.build();
        // All-carbon graphs keep the N/P exclusion heuristics out of play, so
        // the matcher itself is tested against an arbitrary must-set.
        let must_p = 60 + rng.below(41);
        let on = aromatic_atoms(&mol);
        let must: Vec<bool> = (0..n).map(|i| on[i] && rng.below(100) < must_p).collect();
        let result = kekulize_with(&mol, |a| must[a.0 as usize]);
        emit_case(out, &format!("random-{seed}-{case}"), &mol, &must, &[], result)?;
    }
    Ok(())
}

fn smiles_mode(files: &[String], out: &mut impl Write) -> std::io::Result<()> {
    for file in files {
        let bytes = std::fs::read(file)?;
        let text = String::from_utf8_lossy(&bytes);
        for (line_no, line) in text.lines().enumerate() {
            let Some(smi) = line.split_whitespace().next() else { continue };
            let Ok(mut mol) = chematic_smiles::parse(smi) else { continue };
            let mut on = aromatic_atoms(&mol);
            if !on.iter().any(|&a| a) {
                // Kekulé input: perceive aromaticity first, so the
                // aromatize -> kekulize round trip is exercised too.
                chematic_perception::aromatize(&mut mol);
                on = aromatic_atoms(&mol);
                if !on.iter().any(|&a| a) {
                    continue;
                }
            }
            let must: Vec<bool> = (0..mol.atom_count())
                .map(|i| on[i] && atom_must_be_matched(&mol, AtomIdx(i as u32)))
                .collect();
            // Superset of the atoms `kekulize`'s ring-fusion / bridgehead
            // fallbacks may drop from the must-set: N or P with >= 3 aromatic bonds.
            let excl: Vec<u32> = (0..mol.atom_count())
                .filter(|&i| {
                    let idx = AtomIdx(i as u32);
                    must[i]
                        && matches!(mol.atom(idx).element.atomic_number(), 7 | 15)
                        && mol
                            .neighbors(idx)
                            .filter(|&(_, bd)| mol.bond(bd).order == BondOrder::Aromatic)
                            .count()
                            >= 3
                })
                .map(|i| i as u32)
                .collect();
            let result = kekulize(&mol);
            let name = std::path::Path::new(file)
                .file_name()
                .map(|s| s.to_string_lossy().into_owned())
                .unwrap_or_default();
            emit_case(out, &format!("{name}:{}", line_no + 1), &mol, &must, &excl, result)?;
        }
    }
    Ok(())
}

fn main() -> std::io::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let mut out = BufWriter::new(std::io::stdout().lock());
    match args.get(1).map(String::as_str) {
        Some("random") => {
            let seed = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(1);
            let count = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(10_000);
            random_mode(seed, count, &mut out)
        }
        Some("smiles") => smiles_mode(&args[2..], &mut out),
        _ => {
            eprintln!("usage: kekule-diff-harness (random <seed> <count> | smiles <file>...)");
            std::process::exit(2)
        }
    }
}
