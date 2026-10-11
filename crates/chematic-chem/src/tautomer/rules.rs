/// Bond order match type for tautomer rule patterns.
#[derive(Clone, Copy, PartialEq, Eq)]
enum BondOrderMatch {
    Single,
    Double,
    #[allow(dead_code)]
    Any,
}

impl BondOrderMatch {
    fn matches(self, order: BondOrder) -> bool {
        match self {
            BondOrderMatch::Single => {
                matches!(order, BondOrder::Single | BondOrder::Up | BondOrder::Down)
            }
            BondOrderMatch::Double => matches!(order, BondOrder::Double),
            BondOrderMatch::Any => true,
        }
    }
}

/// A tautomer transformation rule: donor loses an H and the bond orders shift.
///
/// For path_len=3 (1,3-shift):
///   Pattern: donor -[donor_bridge_order]- bridge -[bridge_acceptor_order]- acceptor
///   After: donor =[new]= bridge -[new]- acceptor
///
/// For path_len=5 (1,5-shift):
///   Pattern: donor -[donor_bridge_order]- b1 -[any]- b2 -[any]- b3 -[bridge_acceptor_order]- acceptor
///   After: donor =[new]= b1 -[any]- b2 -[any]- b3 -[new]- acceptor (with H shifted)
struct TautomerRule {
    #[allow(dead_code)]
    name: &'static str,
    /// Atomic number of the donor atom (loses H).
    donor_elem: u8,
    /// Atomic number of the bridge atom(s). For 1,3-shift: single bridge. For 1,5-shift: central atom (b2) if specified.
    bridge_elem: Option<u8>,
    /// Atomic number of the acceptor atom (gains H via implicit valence).
    acceptor_elem: u8,
    /// Required bond order between donor and bridge (or donor and b1 for 1,5-shift).
    donor_bridge_order: BondOrderMatch,
    /// Required bond order between bridge and acceptor (or b3 and acceptor for 1,5-shift).
    bridge_acceptor_order: BondOrderMatch,
    /// If true, this rule is applied in canonical_tautomer to normalize toward a preferred form.
    prefer_forward: bool,
    /// Path length: 3 for 1,3-shift (donor-bridge-acceptor), 5 for 1,5-shift.
    path_len: usize,
}

/// The 44 tautomer rules.
static RULES: &[TautomerRule] = &[
    // 1. keto-enol: O-H adjacent to C=C → O=C-C (prefer keto)
    TautomerRule {
        name: "keto-enol",
        donor_elem: 8,
        bridge_elem: Some(6),
        acceptor_elem: 6,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: true,
        path_len: 3,
    },
    // 2. amide-iminol: N-H adjacent to C=O → N=C-O (prefer amide — forward=false)
    TautomerRule {
        name: "amide-iminol",
        donor_elem: 7,
        bridge_elem: Some(6),
        acceptor_elem: 8,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 3,
    },
    // 3. iminol→amide: O-H adjacent to C=N → O=C-N (prefer amide — forward=true)
    TautomerRule {
        name: "iminol-amide",
        donor_elem: 8,
        bridge_elem: Some(6),
        acceptor_elem: 7,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: true,
        path_len: 3,
    },
    // 4. imine-enamine: N-H adjacent to C=C → N=C-C (prefer imine)
    TautomerRule {
        name: "imine-enamine",
        donor_elem: 7,
        bridge_elem: Some(6),
        acceptor_elem: 6,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: true,
        path_len: 3,
    },
    // 5. 1,3-H-shift N→O (any bridge): e.g. nitroso/oxime, hydroxamic acid
    TautomerRule {
        name: "1,3-N-to-O",
        donor_elem: 7,
        bridge_elem: None,
        acceptor_elem: 8,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 3,
    },
    // 6. 1,3-H-shift N→N (any bridge): imidazole, pyrazole, guanidine
    TautomerRule {
        name: "1,3-N-to-N",
        donor_elem: 7,
        bridge_elem: None,
        acceptor_elem: 7,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 3,
    },
    // 7. thioamide: N-H adjacent to C=S → N=C-S (prefer thioamide — forward=false)
    TautomerRule {
        name: "thioamide",
        donor_elem: 7,
        bridge_elem: Some(6),
        acceptor_elem: 16,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 3,
    },
    // 8. thio-iminol→thioamide: S-H adjacent to C=N → S=C-N
    TautomerRule {
        name: "thio-iminol-amide",
        donor_elem: 16,
        bridge_elem: Some(6),
        acceptor_elem: 7,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: true,
        path_len: 3,
    },
    // 9. thio keto-enol: S-H adjacent to C=C → S=C-C (prefer thioketone)
    TautomerRule {
        name: "thio-keto-enol",
        donor_elem: 16,
        bridge_elem: Some(6),
        acceptor_elem: 6,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: true,
        path_len: 3,
    },
    // 10. thio-enol→thioketone: O-H adjacent to C=S → O=C-S
    TautomerRule {
        name: "thio-enol-ketone",
        donor_elem: 8,
        bridge_elem: Some(6),
        acceptor_elem: 16,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 3,
    },
    // 11. 1,3-N→S (any bridge): sulfonamide/thioamide-type
    TautomerRule {
        name: "1,3-N-to-S",
        donor_elem: 7,
        bridge_elem: None,
        acceptor_elem: 16,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 3,
    },
    // 12. 1,3-S→O (any bridge)
    TautomerRule {
        name: "1,3-S-to-O",
        donor_elem: 16,
        bridge_elem: None,
        acceptor_elem: 8,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 3,
    },
    // 13. 1,3-S→N (any bridge)
    TautomerRule {
        name: "1,3-S-to-N",
        donor_elem: 16,
        bridge_elem: None,
        acceptor_elem: 7,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 3,
    },
    // 14. 1,3-O→S (any bridge)
    TautomerRule {
        name: "1,3-O-to-S",
        donor_elem: 8,
        bridge_elem: None,
        acceptor_elem: 16,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 3,
    },
    // 15. 1,3-S→S (any bridge): dithioamide, xanthate-type
    TautomerRule {
        name: "1,3-S-to-S",
        donor_elem: 16,
        bridge_elem: None,
        acceptor_elem: 16,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 3,
    },
    // 16. 1,3-O→N any bridge: extends iminol-amide (rule 3) to non-C bridges
    //     e.g. O-H + S=N, N=N → O=X + N-H (hydroxamic acid, thiohydroximate)
    TautomerRule {
        name: "1,3-O-to-N-any-bridge",
        donor_elem: 8,
        bridge_elem: None,
        acceptor_elem: 7,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: true,
        path_len: 3,
    },
    // 17. 1,3-O→O any bridge: O-H via N or S bridge to =O
    //     e.g. hydroxylamine HO-N=O ↔ O=N-OH (N-oxide tautomer)
    TautomerRule {
        name: "1,3-O-to-O-any-bridge",
        donor_elem: 8,
        bridge_elem: None,
        acceptor_elem: 8,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 3,
    },
    // 18. 1,3-N→C any bridge: extends imine-enamine (rule 4) to non-C bridges
    //     e.g. N-H + S=C, N=C via N bridge → N= + X-C-H
    TautomerRule {
        name: "1,3-N-to-C-any-bridge",
        donor_elem: 7,
        bridge_elem: None,
        acceptor_elem: 6,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: true,
        path_len: 3,
    },
    // 19. 1,3-C→O any bridge: active methylene (or terminal alkyne) to O-H
    //     Forward: C-H + X=O → C=X + O-H. prefer_forward:false → prefer keto/C-H form.
    TautomerRule {
        name: "1,3-C-to-O-any-bridge",
        donor_elem: 6,
        bridge_elem: None,
        acceptor_elem: 8,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 3,
    },
    // 20. nitroso → oxime: C-H adjacent to N=O → C=N-OH.
    // Keep this narrower than the generic C→O rule above: enabling that
    // rule in the forward direction would incorrectly enolize ordinary
    // carbonyls through arbitrary bridges.  This is the concrete
    // nitroso/oxime case documented in Phase 2 section 1.6.
    TautomerRule {
        name: "nitroso-oxime",
        donor_elem: 6,
        bridge_elem: Some(7),
        acceptor_elem: 8,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: true,
        path_len: 3,
    },
    // 21. 1,3-C→N any bridge: active methylene adjacent to =N (via S, O, or N bridge)
    //     Forward: C-H + X=N → C=X + N-H. prefer_forward:false → prefer N-H form.
    TautomerRule {
        name: "1,3-C-to-N-any-bridge",
        donor_elem: 6,
        bridge_elem: None,
        acceptor_elem: 7,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 3,
    },
    // 1,5-H shift rules (path_len=5)
    // 22. 1,5-O→O: β-diketone (acetylacetone)
    //     Pattern: O-C(-)-C(-)-C(=O) → O=C(-)-C(-)-C(-O)
    TautomerRule {
        name: "1,5-O-to-O-beta-diketone",
        donor_elem: 8,
        bridge_elem: Some(6), // Central carbon in pattern (b2)
        acceptor_elem: 8,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 5,
    },
    // 23. 1,5-O→N: enol imine
    TautomerRule {
        name: "1,5-O-to-N",
        donor_elem: 8,
        bridge_elem: Some(6),
        acceptor_elem: 7,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 5,
    },
    // 24. 1,5-N→N: extended guanidine/amidine tautomerism
    TautomerRule {
        name: "1,5-N-to-N",
        donor_elem: 7,
        bridge_elem: Some(6),
        acceptor_elem: 7,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 5,
    },
    // 25. 1,5-N→O: hydroxamic acid type
    TautomerRule {
        name: "1,5-N-to-O",
        donor_elem: 7,
        bridge_elem: Some(6),
        acceptor_elem: 8,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 5,
    },
    // 26. 1,5-C→O: active methylene with conjugation
    TautomerRule {
        name: "1,5-C-to-O",
        donor_elem: 6,
        bridge_elem: Some(6),
        acceptor_elem: 8,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 5,
    },
    // 27. 1,5-O→O with N bridge: nitro-type tautomerism
    //     e.g. O-C-N(=O)-C-O ↔ O=C-N(-O)-C-O via N bridge
    TautomerRule {
        name: "1,5-O-to-O-N-bridge",
        donor_elem: 8,
        bridge_elem: Some(7),
        acceptor_elem: 8,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 5,
    },
    // 28. 1,5-O→O with S bridge: thio-β-diketone
    TautomerRule {
        name: "1,5-O-to-O-S-bridge",
        donor_elem: 8,
        bridge_elem: Some(16),
        acceptor_elem: 8,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 5,
    },
    // 29. 1,5-N→O with C bridge (existing, carbon specified)
    // Already covered by rule 22 (path_len=5 with C bridge)
    // 30. 1,5-N→O with N bridge: bridging N (amidino-type)
    TautomerRule {
        name: "1,5-N-to-O-N-bridge",
        donor_elem: 7,
        bridge_elem: Some(7),
        acceptor_elem: 8,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 5,
    },
    // 31. 1,5-N→O with S bridge
    TautomerRule {
        name: "1,5-N-to-O-S-bridge",
        donor_elem: 7,
        bridge_elem: Some(16),
        acceptor_elem: 8,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 5,
    },
    // 32. 1,5-N→N with N bridge: guanidine-type via N
    TautomerRule {
        name: "1,5-N-to-N-N-bridge",
        donor_elem: 7,
        bridge_elem: Some(7),
        acceptor_elem: 7,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 5,
    },
    // 33. 1,5-N→N with S bridge
    TautomerRule {
        name: "1,5-N-to-N-S-bridge",
        donor_elem: 7,
        bridge_elem: Some(16),
        acceptor_elem: 7,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 5,
    },
    // 34. 1,5-O→N with N bridge: hydroxamic-type via N
    TautomerRule {
        name: "1,5-O-to-N-N-bridge",
        donor_elem: 8,
        bridge_elem: Some(7),
        acceptor_elem: 7,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 5,
    },
    // 35. 1,5-O→N with S bridge
    TautomerRule {
        name: "1,5-O-to-N-S-bridge",
        donor_elem: 8,
        bridge_elem: Some(16),
        acceptor_elem: 7,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 5,
    },
    // 36. 1,5-S→O with N bridge
    TautomerRule {
        name: "1,5-S-to-O-N-bridge",
        donor_elem: 16,
        bridge_elem: Some(7),
        acceptor_elem: 8,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 5,
    },
    // 37. 1,5-S→O with S bridge
    TautomerRule {
        name: "1,5-S-to-O-S-bridge",
        donor_elem: 16,
        bridge_elem: Some(16),
        acceptor_elem: 8,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 5,
    },
    // 38. 1,5-S→N with C bridge
    TautomerRule {
        name: "1,5-S-to-N-C-bridge",
        donor_elem: 16,
        bridge_elem: Some(6),
        acceptor_elem: 7,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 5,
    },
    // 39. 1,5-S→N with N bridge
    TautomerRule {
        name: "1,5-S-to-N-N-bridge",
        donor_elem: 16,
        bridge_elem: Some(7),
        acceptor_elem: 7,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 5,
    },
    // 40. 1,5-C→N with C bridge: extended enamine-imine
    TautomerRule {
        name: "1,5-C-to-N-C-bridge",
        donor_elem: 6,
        bridge_elem: Some(6),
        acceptor_elem: 7,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 5,
    },
    // 41. 1,5-C→N with N bridge
    TautomerRule {
        name: "1,5-C-to-N-N-bridge",
        donor_elem: 6,
        bridge_elem: Some(7),
        acceptor_elem: 7,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 5,
    },
    // 42. 1,5-C→N with S bridge
    TautomerRule {
        name: "1,5-C-to-N-S-bridge",
        donor_elem: 6,
        bridge_elem: Some(16),
        acceptor_elem: 7,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 5,
    },
    // 43. 1,5-C→S with N bridge
    TautomerRule {
        name: "1,5-C-to-S-N-bridge",
        donor_elem: 6,
        bridge_elem: Some(7),
        acceptor_elem: 16,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 5,
    },
    // 44. 1,5-C→S with S bridge
    TautomerRule {
        name: "1,5-C-to-S-S-bridge",
        donor_elem: 6,
        bridge_elem: Some(16),
        acceptor_elem: 16,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 5,
    },
];
