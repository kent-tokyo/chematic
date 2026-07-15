//! Milestone 4B-1.5: the pure-operator metamorphic suite the user specified, run
//! *before* trusting any molecule-level result. These test the Like/Unlike comparator
//! itself on synthetic `DescriptorFamily` sequences -- no digraph, no molecule,
//! Like-wins pinned exactly as the user's spec states ("最初に現れるlike descriptor
//! pairが、対応するunlike pairに優先する"). If a molecule-level gate needs a different
//! winner polarity to reach 8/8, the bug is elsewhere (reference selection, or the
//! rank/label mapping), not in this operator -- these tests exist to make that
//! distinction possible instead of guessed.
//!
//! Usage: cargo run -p chematic-cip --release --example rule4b_operator_tests

use std::cmp::Ordering;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DescriptorFamily {
    Right, // R, M, seqCis
    Left,  // S, P, seqTrans
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PairRelation {
    Like,
    Unlike,
}

fn relation(reference: DescriptorFamily, value: DescriptorFamily) -> PairRelation {
    if value == reference {
        PairRelation::Like
    } else {
        PairRelation::Unlike
    }
}

/// Compare two sequences against the SAME reference. `Greater` means `a` outranks `b`.
/// Like precedes Unlike at the first differing position -- pinned per the user's spec,
/// never varied elsewhere in this file.
fn compare_same_reference(
    reference: DescriptorFamily,
    a: &[DescriptorFamily],
    b: &[DescriptorFamily],
) -> Ordering {
    for (&fa, &fb) in a.iter().zip(b.iter()) {
        let (ra, rb) = (relation(reference, fa), relation(reference, fb));
        match (ra, rb) {
            (PairRelation::Like, PairRelation::Unlike) => return Ordering::Greater,
            (PairRelation::Unlike, PairRelation::Like) => return Ordering::Less,
            _ => {}
        }
    }
    Ordering::Equal
}

/// Compare two sequences against their OWN, independently-chosen references -- the
/// user's actual spec (never a shared cross-branch reference). Still Like-wins.
fn compare_independent_references(
    ref_a: DescriptorFamily,
    seq_a: &[DescriptorFamily],
    ref_b: DescriptorFamily,
    seq_b: &[DescriptorFamily],
) -> Ordering {
    for (&fa, &fb) in seq_a.iter().zip(seq_b.iter()) {
        let (ra, rb) = (relation(ref_a, fa), relation(ref_b, fb));
        match (ra, rb) {
            (PairRelation::Like, PairRelation::Unlike) => return Ordering::Greater,
            (PairRelation::Unlike, PairRelation::Like) => return Ordering::Less,
            _ => {}
        }
    }
    Ordering::Equal
}

/// Candidate-count rule: fewer reference candidates outranks more (a ligand with an
/// unambiguous single reference beats one with an ambiguous tied pair) -- `None` if
/// counts are equal (falls through to comparing the candidate pair lists themselves,
/// not modeled by this simple helper).
fn compare_candidate_counts(a_count: usize, b_count: usize) -> Option<Ordering> {
    match a_count.cmp(&b_count) {
        Ordering::Equal => None,
        Ordering::Less => Some(Ordering::Greater),
        Ordering::Greater => Some(Ordering::Less),
    }
}

/// Reference-candidate selection: majority family among descriptors at the first
/// level; both families if the count is tied.
fn reference_candidates(descriptors: &[DescriptorFamily]) -> Vec<DescriptorFamily> {
    let right = descriptors
        .iter()
        .filter(|&&f| f == DescriptorFamily::Right)
        .count();
    let left = descriptors.len() - right;
    match right.cmp(&left) {
        Ordering::Greater => vec![DescriptorFamily::Right],
        Ordering::Less => vec![DescriptorFamily::Left],
        Ordering::Equal => vec![DescriptorFamily::Right, DescriptorFamily::Left],
    }
}

struct TestResult {
    name: &'static str,
    passed: bool,
    detail: String,
}

fn check(name: &'static str, passed: bool, detail: impl Into<String>) -> TestResult {
    TestResult {
        name,
        passed,
        detail: detail.into(),
    }
}

fn main() {
    use DescriptorFamily::{Left, Right};

    let mut results = Vec::new();

    // 1. RR > RS (reference = R).
    let ord = compare_same_reference(Right, &[Right, Right], &[Right, Left]);
    results.push(check(
        "RR > RS (ref=R)",
        ord == Ordering::Greater,
        format!("compare_same_reference(R, [R,R], [R,S]) = {ord:?}, want Greater"),
    ));

    // 2. SS > SR (reference = S).
    let ord = compare_same_reference(Left, &[Left, Left], &[Left, Right]);
    results.push(check(
        "SS > SR (ref=S)",
        ord == Ordering::Greater,
        format!("compare_same_reference(S, [S,S], [S,R]) = {ord:?}, want Greater"),
    ));

    // 3. Mirror (global inversion) invariance: flipping reference AND every value
    // together must preserve which sequence wins.
    let ord_original = compare_same_reference(Right, &[Right, Right], &[Right, Left]);
    let ord_mirrored = compare_same_reference(Left, &[Left, Left], &[Left, Right]);
    results.push(check(
        "mirror invariance (RR>RS mirrors to SS>SR, same winner)",
        ord_original == ord_mirrored,
        format!("original={ord_original:?} mirrored={ord_mirrored:?}"),
    ));

    // 4. Branch antisymmetry: compare(a,b) == compare(b,a).reverse().
    let forward = compare_same_reference(Right, &[Right, Left], &[Left, Right]);
    let backward = compare_same_reference(Right, &[Left, Right], &[Right, Left]);
    results.push(check(
        "branch antisymmetry: compare(a,b) == reverse(compare(b,a))",
        forward == backward.reverse(),
        format!("forward={forward:?} backward={backward:?}"),
    ));

    // 5. Independent references (left ref=S, right ref=R -- different families):
    // left=[S,S] (all Like relative to its own ref S) vs right=[R,S] (Like,Unlike
    // relative to its own ref R) -- left should win, decided at position 1.
    let ord = compare_independent_references(Left, &[Left, Left], Right, &[Right, Left]);
    results.push(check(
        "independent references: left(ref S)=[S,S] > right(ref R)=[R,S]",
        ord == Ordering::Greater,
        format!("compare_independent_references(...) = {ord:?}, want Greater"),
    ));

    // 6. Candidate-count rule: 1 candidate outranks 2.
    let ord = compare_candidate_counts(1, 2);
    results.push(check(
        "candidate-count: 1 candidate > 2 candidates",
        ord == Some(Ordering::Greater),
        format!("compare_candidate_counts(1,2) = {ord:?}, want Some(Greater)"),
    ));

    // 7. Equal-majority: R/S tied count produces BOTH candidates, order-independent
    // (checked both input orders give the same *set*, since Vec order from a stable
    // count-based construction shouldn't matter for this property -- comparing as
    // unordered).
    let mut a = reference_candidates(&[Right, Left]);
    let mut b = reference_candidates(&[Left, Right]);
    a.sort_by_key(|f| matches!(f, Right) as u8);
    b.sort_by_key(|f| matches!(f, Right) as u8);
    results.push(check(
        "equal-majority (1R/1S) produces both candidates, order-independent",
        a == vec![Left, Right] && a == b,
        format!("[R,S] input -> {a:?}; [S,R] input -> {b:?}"),
    ));

    let mut all_pass = true;
    for r in &results {
        let mark = if r.passed { "PASS" } else { "FAIL" };
        println!("[{mark}] {}: {}", r.name, r.detail);
        all_pass &= r.passed;
    }
    println!(
        "\n{}/{} operator tests passed.",
        results.iter().filter(|r| r.passed).count(),
        results.len()
    );
    if !all_pass {
        std::process::exit(1);
    }

    println!(
        "\nNote: 'reference-dependent reordering' (an artificial digraph whose child \
         priority order changes depending on which reference is fixed) is NOT tested \
         here -- it requires the Step-B re-sorter, not yet built. Left as an open item."
    );
}
