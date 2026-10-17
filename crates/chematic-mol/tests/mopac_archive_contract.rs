use chematic_mol::{MopacParseLimits, parse_mopac_archive, parse_mopac_archive_with_limits};

const INTERNAL: &str = "Summary\nFINAL GEOMETRY OBTAINED\nPM7 SYMMETRY\nSynthetic formaldehyde\n\nO 0 +0 0 +0 0 +0 0 0 0\nC 1.2 +1 0 +0 0 +0 1 0 0\nH 1.1 +1 120 +1 0 +0 2 1 0\nH 1.1 +0 120 +0 180 +0 2 1 3\n\n3 1 4\n";
const CART: &str = "Summary\n FINAL GEOMETRY OBTAINED\nPM7 XYZ\nWater\nCartesian\nO 0 +1 0 +1 0 +1\nH 1 +1 0 +1 0 +1\nH 0 +1 1 +1 0 +1\n";

#[test]
fn internal_values_references_and_ancillary_sections_are_preserved() {
    let a = parse_mopac_archive(INTERNAL).unwrap();
    assert_eq!(a.atoms.len(), 4);
    assert_eq!(a.atoms[3].values, [1.1, 120., 180.]);
    assert_eq!(a.atoms[3].references, Some([2, 1, 3]));
    assert_eq!(a.summary_text, "Summary\n");
    assert_eq!(a.trailing_text, "\n3 1 4\n");
    assert_eq!(a.source_text, INTERNAL);
    assert_eq!(a.restart_text, INTERNAL.split_once("OBTAINED\n").unwrap().1);
}

#[test]
fn cartesian_values_are_not_labeled_internal() {
    let a = parse_mopac_archive(CART).unwrap();
    assert!(a.atoms.iter().all(|a| a.references.is_none()));
    assert_eq!(a.atoms[2].values, [0., 1., 0.]);
    assert_eq!(a.trailing_text, "");
}

#[test]
fn charge_column_requires_announcement_and_is_finite() {
    let text = CART
        .replace("OBTAINED", "OBTAINED CHARGE")
        .replace("0 +1\n", "0 +1 -0.25\n");
    let a = parse_mopac_archive(&text).unwrap();
    assert!(a.atoms.iter().all(|a| a.partial_charge == Some(-0.25)));
    assert!(parse_mopac_archive(&text.replace("OBTAINED CHARGE", "OBTAINED")).is_err());
    assert!(parse_mopac_archive(&text.replace("-0.25", "NaN")).is_err());
}

#[test]
fn supports_crlf_fortran_exponents_labels_and_atomic_numbers() {
    let text = CART.replace("H 1", "1(edge) 1D0").replace('\n', "\r\n");
    let a = parse_mopac_archive(&text).unwrap();
    assert_eq!(a.atoms[1].element.symbol(), "H");
    assert_eq!(a.atoms[1].label.as_deref(), Some("edge"));
    assert_eq!(a.source_text, text);
}

#[test]
fn malformed_or_ambiguous_archives_fail_without_partial_results() {
    for text in [
        CART.replace("FINAL GEOMETRY OBTAINED", ""),
        format!("{CART}{CART}"),
        CART.replace("OBTAINED", "OBTAINED UNKNOWN"),
        CART.replace("O 0 +1 0 +1 0 +1", ""),
        CART.replace("O 0 +1", "O NaN +1"),
        CART.replace("O 0 +1", "O 0 3"),
        CART.replace("PM7 XYZ", "A0"),
        CART.replace("PM7 XYZ", "PM7 +"),
        CART.replace("O 0 +1", "Tv 0 +1"),
        CART.replace("O 0 +1", "99 0 +1"),
        INTERNAL.replace("2 1 3", "2 1 4"),
        INTERNAL.replace("2 1 3", "2 1 5"),
        INTERNAL.replace("2 1 3", "2 1 2"),
        INTERNAL.replace("2 1 3", "2 1 -3"),
        CART.replace("Summary", "Sum\0mary"),
    ] {
        assert!(parse_mopac_archive(&text).is_err(), "{text}");
    }
}

#[test]
fn archive_error_uses_original_line_number() {
    let err = parse_mopac_archive(&CART.replace("H 1 +1", "H NaN +1")).unwrap_err();
    assert_eq!(err.line, 7);
}

#[test]
fn all_limits_are_enforced_on_the_whole_archive() {
    for limits in [
        MopacParseLimits {
            max_input_bytes: CART.len() - 1,
            ..Default::default()
        },
        MopacParseLimits {
            max_atoms: 2,
            ..Default::default()
        },
        MopacParseLimits {
            max_line_len: 10,
            ..Default::default()
        },
    ] {
        assert!(parse_mopac_archive_with_limits(CART, &limits).is_err());
    }
    assert!(
        parse_mopac_archive_with_limits(
            CART,
            &MopacParseLimits {
                max_input_bytes: CART.len(),
                max_atoms: 3,
                max_line_len: 24,
            }
        )
        .is_ok()
    );
}
