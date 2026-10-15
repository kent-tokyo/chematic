from pathlib import Path
r=Path('/private/tmp/chematic-openbabel-ring-queue-20261015');s=Path('/private/tmp/chematic-openbabel-ring-queue-source-20261015');p=s/'crates/chematic-perception/src/stereo2d_ez_direction.rs';old=p.read_text();a=old.index('fn in_ring_smaller_than(');b=old.index('fn has_other_double_bond(',a);reference=old[a:b].replace('fn in_ring_smaller_than(', 'fn reference_ring_search(',1).replace('-> bool {','-> (bool, usize) {',1).replace('return true;','return (true, seen.len());',1);i=reference.rfind('    false');reference=reference[:i]+reference[i:].replace('    false','    (false, seen.len())',1)
new='''fn in_ring_smaller_than(mol: &Molecule, bond: BondIdx, size: usize) -> bool {
    const INLINE: usize = 32;
    let (start, goal) = (mol.bond(bond).atom1, mol.bond(bond).atom2);
    // The queue is also the complete visited list. Small neighborhoods need
    // no heap allocation; larger ones keep hashed membership before spilling.
    let mut queue: smallvec::SmallVec<[AtomIdx; INLINE]> = smallvec::smallvec![start];
    let mut overflow_seen: Option<HashSet<AtomIdx>> = None;
    let mut begin = 0;
    for _ in 0..size.saturating_sub(2) {
        let end = queue.len();
        for i in begin..end {
            let a = queue[i];
            for (nb, bi) in mol.neighbors(a) {
                if bi == bond {
                    continue;
                }
                if nb == goal {
                    return true;
                }
                let fresh = match &mut overflow_seen {
                    Some(seen) => seen.insert(nb),
                    None => {
                        if queue.contains(&nb) {
                            false
                        } else {
                            if queue.len() == INLINE {
                                let mut seen: HashSet<_> = queue.iter().copied().collect();
                                seen.insert(nb);
                                overflow_seen = Some(seen);
                            }
                            true
                        }
                    }
                };
                if fresh {
                    queue.push(nb);
                }
            }
        }
        begin = end;
        if begin == queue.len() {
            break;
        }
    }
    false
}

'''
p.write_text(old[:a]+new+old[b:]+'\n#[cfg(test)]\n#[path = "stereo2d_ez_ring_queue_tests.rs"]\nmod ring_queue_tests;\n')
test='//! Baseline breadth-first traversal with return-only visited-count instrumentation.\nuse super::*;\n\n'+reference+'''
#[test]
fn ring_queue_keeps_original_depth_and_membership_on_all_pinned_bonds() {
    let corpora = [include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi")),include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../validation/cosmolkit_comparison/corpora/stress-v1.smi"))];
    let sizes=(0..=12).chain([33,65]).collect::<Vec<_>>();
    let mut rows=0usize;
    let mut queries=0usize;
    let mut inline_cases=0usize;
    let mut hashed_cases=0usize;
    let mut found=0usize;
    let mut not_found=0usize;
    for corpus in corpora {
        for line in corpus.lines().filter(|line| !line.trim().is_empty()) {
            let mol=chematic_smiles::parse_template(line.split_whitespace().next().unwrap()).unwrap();
            for (bi,_) in mol.bonds() {
                for &size in &sizes {
                    let (expected,visited)=reference_ring_search(&mol,bi,size);
                    assert_eq!(in_ring_smaller_than(&mol,bi,size),expected,"row {rows}, bond {bi:?}, size {size}");
                    queries+=1;
                    if visited<=32 { inline_cases+=1; } else { hashed_cases+=1; }
                    if expected { found+=1; } else { not_found+=1; }
                }
            }
            rows+=1;
        }
    }
    assert_eq!(rows,14072);
    assert!(queries>1_000_000);
    assert!(inline_cases>0 && hashed_cases>0 && found>0 && not_found>0);
    eprintln!("ring_queue_oracle {}",serde_json::json!({"rows":rows,"queries":queries,"inline_cases":inline_cases,"hashed_cases":hashed_cases,"found":found,"not_found":not_found}));
}
''';(s/'crates/chematic-perception/src/stereo2d_ez_ring_queue_tests.rs').write_text(test)
p=s/'CHANGELOG.md';p.write_text(p.read_text().replace('### Performance\n','### Performance\n\n- Keep E/Z small-ring searches in one inline visited queue for up to 32 atoms;\n  retain hashed membership for larger neighborhoods and the same breadth-first\n  depth boundary, avoiding per-layer temporary arrays.\n',1));print('Inline ring-search queue with hashed overflow; exact baseline BFS oracle written')
