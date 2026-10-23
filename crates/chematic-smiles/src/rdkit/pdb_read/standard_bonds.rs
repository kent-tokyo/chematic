//! Ordered PDB residue bond-order rules from `StandardPDBDoubleBond`.

fn atom_pair_is(a1: [u8; 4], a2: [u8; 4], x: &[u8; 4], y: &[u8; 4]) -> bool {
    a1 == *x && a2 == *y
}

fn is_backbone_residue(res: &[u8]) -> bool {
    matches!(
        res,
        b"ACE"
            | b"ALA"
            | b"CYS"
            | b"GLY"
            | b"ILE"
            | b"LEU"
            | b"LYS"
            | b"MET"
            | b"PRO"
            | b"SER"
            | b"THR"
            | b"VAL"
    )
}

fn is_protein_residue(res: &[u8]) -> bool {
    is_backbone_residue(res)
        || matches!(
            res,
            b"ARG" | b"ASN" | b"ASP" | b"GLN" | b"GLU" | b"HIS" | b"PHE" | b"TYR" | b"TRP"
        )
}

fn protein_sidechain_double_bond(res: &[u8], a1: [u8; 4], a2: [u8; 4]) -> bool {
    let is = |x, y| atom_pair_is(a1, a2, x, y);
    match res {
        b"ARG" => is(b" CZ ", b" NH2"),
        b"ASN" | b"ASP" => is(b" CG ", b" OD1"),
        b"GLN" | b"GLU" => is(b" CD ", b" OE1"),
        b"HIS" => is(b" CD2", b" CG ") || is(b" CE1", b" ND1"),
        b"PHE" | b"TYR" => is(b" CD1", b" CG ") || is(b" CD2", b" CE2") || is(b" CE1", b" CZ "),
        b"TRP" => {
            is(b" CD1", b" CG ")
                || is(b" CD2", b" CE2")
                || is(b" CE3", b" CZ3")
                || is(b" CH2", b" CZ2")
        }
        _ => false,
    }
}

fn purine_double_bond(a1: [u8; 4], a2: [u8; 4], guanine: bool) -> bool {
    let is = |x, y| atom_pair_is(a1, a2, x, y);
    (if guanine {
        is(b" C6 ", b" O6 ")
    } else {
        is(b" C6 ", b" N1 ")
    }) || is(b" C2 ", b" N3 ")
        || is(b" C4 ", b" C5 ")
        || is(b" C8 ", b" N7 ")
        || is(b" OP1", b" P  ")
}

fn pyrimidine_double_bond(a1: [u8; 4], a2: [u8; 4], cytosine: bool) -> bool {
    let is = |x, y| atom_pair_is(a1, a2, x, y);
    is(b" C2 ", b" O2 ")
        || if cytosine {
            is(b" C4 ", b" N3 ") || is(b" C5 ", b" C6 ") || is(b" OP1", b" P  ")
        } else {
            is(b" C4 ", b" O4 ") || is(b" C5 ", b" C6 ") || is(b" OP1", b" P  ")
        }
}

fn nucleic_double_bond(res: &[u8], a1: [u8; 4], a2: [u8; 4]) -> bool {
    match res {
        b"  A" | b" DA" => purine_double_bond(a1, a2, false),
        b"  G" | b" DG" => purine_double_bond(a1, a2, true),
        b"  C" | b" DC" => pyrimidine_double_bond(a1, a2, true),
        b"  T" | b" DT" | b"  U" | b" DU" => pyrimidine_double_bond(a1, a2, false),
        _ => false,
    }
}

/// Whether a standard residue atom pair is double bonded.
pub(super) fn is_double_bond(res: &[u8], a1: [u8; 4], a2: [u8; 4]) -> bool {
    let (a1, a2) = if u32::from_be_bytes(a1) > u32::from_be_bytes(a2) {
        (a2, a1)
    } else {
        (a1, a2)
    };
    let backbone = atom_pair_is(a1, a2, b" C  ", b" O  ");
    (backbone && is_protein_residue(res))
        || protein_sidechain_double_bond(res, a1, a2)
        || nucleic_double_bond(res, a1, a2)
}
