//! RDKit's valence model (RDKit 2026.03), used where reaction products must
//! agree with what RDKit's `SanitizeMol` accepts and with the implicit
//! hydrogen counts it assigns (issue #734).
//!
//! RDKit looks a charged atom up as its isoelectronic element (`N+` as `C`,
//! `O-` as `F`, `Cl+` as `S`), except that P/S (and As/Se) pushed past the
//! next noble-gas-like element keep their own, hypervalent list shifted by the
//! charge (`S-` allows 5, `P-2` allows 3). Elements whose list RDKit leaves
//! open (alkali and alkaline-earth metals, transition metals, lanthanides…)
//! accept any valence.

/// RDKit's allowed valences for element `z` (ascending); `None` = any.
fn valence_list(z: u8) -> Option<&'static [u8]> {
    Some(match z {
        1 => &[1],
        2 | 10 | 18 | 36 | 86 => &[0],
        4 => &[2],
        5 | 13 | 31 | 49 => &[3],
        6 | 14 | 32 => &[4],
        7 => &[3],
        8 => &[2],
        9 | 17 | 35 | 55 | 87 => &[1],
        15 | 33 | 51 | 83 => &[3, 5],
        16 | 34 | 52 | 84 => &[2, 4, 6],
        50 | 82 => &[2, 4],
        53 | 85 => &[1, 3, 5],
        54 => &[0, 2, 4, 6],
        _ => return None,
    })
}

/// The valence list RDKit applies to element `z` with formal `charge`, as
/// `(list, offset)`: each allowed valence is `list[i] + offset`. `None` means
/// RDKit accepts any valence for this atom.
pub(crate) fn allowed_valences(z: u8, charge: i8) -> Option<(&'static [u8], i16)> {
    let own = valence_list(z)?;
    let effective = i16::from(z) - i16::from(charge);
    if effective <= 0 {
        return None;
    }
    let effective = effective.min(118) as u8;
    // P/S and As/Se past the next "noble gas" keep their own list.
    if (effective > 16 && matches!(z, 15 | 16)) || (effective > 34 && matches!(z, 33 | 34)) {
        return Some((own, i16::from(charge)));
    }
    valence_list(effective).map(|list| (list, 0))
}

/// RDKit's largest allowed explicit valence; `None` = unrestricted.
pub(crate) fn max_valence(z: u8, charge: i8) -> Option<i16> {
    if charge == 0 {
        return valence_list(z).map(|list| i16::from(*list.last().expect("non-empty list")));
    }
    if z == 1 && charge == -1 {
        // Hydride may bridge two atoms.
        return Some(2);
    }
    let (list, offset) = allowed_valences(z, charge)?;
    Some(i16::from(*list.last().expect("non-empty valence list")) + offset)
}

/// RDKit's implicit hydrogen count for an atom of element `z`, formal
/// `charge` and explicit valence `explicit` (bond orders plus explicit H):
/// the gap to the smallest allowed valence at or above `explicit`; zero for
/// unrestricted elements or an over-valent atom.
pub(crate) fn implicit_hydrogens(z: u8, charge: i8, explicit: i16) -> u8 {
    if z == 1 && charge == -1 {
        return 0;
    }
    let Some((list, offset)) = allowed_valences(z, charge) else {
        return 0;
    };
    list.iter()
        .map(|&v| i16::from(v) + offset)
        .find(|&v| v >= explicit)
        .map_or(0, |v| (v - explicit).clamp(0, i16::from(u8::MAX)) as u8)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_rdkit_maximum_valences() {
        // (z, charge, RDKit 2026.03.6 largest accepted H count on a bare atom;
        // None = any), from SanitizeMol on `[XHn±c]`.
        for (z, charge, expected) in [
            (6, 0, Some(4)),
            (6, 1, Some(3)),
            (6, -1, Some(3)),
            (6, 3, None),
            (7, 1, Some(4)),
            (7, 2, Some(3)),
            (7, -3, Some(0)),
            (8, -3, None),
            (5, -1, Some(4)),
            (5, 3, Some(0)),
            (15, 0, Some(5)),
            (15, -1, Some(6)),
            (15, -2, Some(3)),
            (15, -3, Some(2)),
            (16, -1, Some(5)),
            (16, -3, Some(3)),
            (17, 1, Some(6)),
            (17, 3, Some(4)),
            (17, -2, None),
            (53, -1, Some(6)),
            (53, -2, Some(1)),
            (52, -3, Some(1)),
            (13, -3, Some(6)),
            (3, 0, None),
            (12, 3, None),
            (26, 0, None),
            (1, 0, Some(1)),
            (1, -1, Some(2)),
            (1, 1, None),
        ] {
            assert_eq!(max_valence(z, charge), expected, "z={z} charge={charge}");
        }
    }

    #[test]
    fn implicit_hydrogens_follow_rdkit() {
        assert_eq!(implicit_hydrogens(34, 0, 1), 1); // C[SeH]
        assert_eq!(implicit_hydrogens(34, 0, 3), 1); // C[SeH](C)C
        assert_eq!(implicit_hydrogens(14, 0, 1), 3); // C[SiH3]
        assert_eq!(implicit_hydrogens(50, 0, 1), 1); // C[SnH]
        assert_eq!(implicit_hydrogens(34, 1, 1), 2); // C[SeH2+]
        assert_eq!(implicit_hydrogens(17, 1, 3), 1); // C[ClH+](C)C
        assert_eq!(implicit_hydrogens(53, -1, 5), 1); // C[IH-](C)(C)(C)C
        assert_eq!(implicit_hydrogens(78, 0, 1), 0); // C[Pt]
        assert_eq!(implicit_hydrogens(7, 0, 4), 0); // over-valent
    }
}
