//! Invalid input and resource contracts of the approximate pure-Rust reader.
use crate::parser::{InchiParseError, InchiParseLimits, parse_inchi, parse_inchi_with_limits};

#[test]
fn malformed_layers_return_typed_errors_without_partial_graphs() {
    for (source, expected) in [
        ("InChI=1S", InchiParseError::InvalidFormat),
        ("InChI=1S/", InchiParseError::InvalidFormula),
        ("InChI=1S/ch4", InchiParseError::InvalidFormula),
        ("InChI=1S/Xx", InchiParseError::InvalidFormula),
        (
            "InChI=1S/C999999999999999999999999999999",
            InchiParseError::InvalidFormula,
        ),
        (
            "InChI=1S/C2H6/c1-3/h1-2H3",
            InchiParseError::InvalidConnectivity,
        ),
        (
            "InChI=1S/C2H6/c1(3)2/h1-2H3",
            InchiParseError::InvalidConnectivity,
        ),
        ("InChI=1S/CH4/h1", InchiParseError::InvalidHydrogen),
        ("InChI=1S/CH4/h1H256", InchiParseError::InvalidHydrogen),
        ("InChI=1S/CH4/hxH4", InchiParseError::InvalidHydrogen),
        ("InChI=1S/CH4/hx-1H4", InchiParseError::InvalidHydrogen),
        ("InChI=1S/CH4/h1-xH4", InchiParseError::InvalidHydrogen),
        (
            "InChI=1S/CH4/tx+",
            InchiParseError::Unsupported("invalid atom number in stereo layer".into()),
        ),
        (
            "InChI=1S/CH4/tx-",
            InchiParseError::Unsupported("invalid atom number in stereo layer".into()),
        ),
    ] {
        let error = parse_inchi(source)
            .err()
            .unwrap_or_else(|| panic!("accepted {source}"));
        assert_eq!(error, expected, "{source}");
        assert!(!error.to_string().is_empty());
    }
}

#[test]
fn parser_budgets_allow_exact_sizes_and_reject_atom_count_overflow() {
    let source = "InChI=1S/C2H6/c1-2/h1-2H3";
    let exact = InchiParseLimits {
        max_input_bytes: source.len(),
        max_atoms: 2,
    };
    assert_eq!(
        parse_inchi_with_limits(source, &exact)
            .unwrap()
            .atom_count(),
        2
    );
    for limits in [
        InchiParseLimits {
            max_input_bytes: source.len() - 1,
            ..exact
        },
        InchiParseLimits {
            max_atoms: 1,
            ..exact
        },
    ] {
        let error = parse_inchi_with_limits(source, &limits).err().unwrap();
        let InchiParseError::ResourceLimit {
            resource,
            actual,
            limit,
        } = error
        else {
            panic!("{error}")
        };
        assert!(actual > limit);
        assert!(error.to_string().contains(resource));
    }
    let source = format!("InChI=1S/C{}N", usize::MAX);
    assert!(matches!(
        parse_inchi(&source),
        Err(InchiParseError::ResourceLimit {
            resource: "atoms",
            actual: usize::MAX,
            ..
        })
    ));
}
