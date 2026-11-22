use chematic_mol::{parse_dlpoly_config, write_dlpoly_config};

#[test]
fn edited_config_rejects_invalid_titles_cells_and_overlong_headers() {
    let original = parse_dlpoly_config("title\n0 0\nX 1\n0 0 0\n").unwrap();
    let mut invalid = original.clone();
    invalid.title = "水".into();
    assert!(
        write_dlpoly_config(&invalid)
            .unwrap_err()
            .detail
            .contains("title must be ASCII")
    );
    invalid = original.clone();
    invalid.title = "x".repeat(73);
    assert!(
        write_dlpoly_config(&invalid)
            .unwrap_err()
            .detail
            .contains("title must be ASCII")
    );
    invalid = original.clone();
    invalid.levcfg = 3;
    assert!(
        write_dlpoly_config(&invalid)
            .unwrap_err()
            .detail
            .contains("supported levcfg")
    );
    invalid = original.clone();
    invalid.imcon = 1;
    assert!(
        write_dlpoly_config(&invalid)
            .unwrap_err()
            .detail
            .contains("cell vectors must be present")
    );
    invalid.cell_vectors = Some([[0.0; 3], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
    assert!(
        write_dlpoly_config(&invalid)
            .unwrap_err()
            .detail
            .contains("zero-length cell vector")
    );
    invalid.cell_vectors = Some([[1.0, 0.1, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
    assert!(
        write_dlpoly_config(&invalid)
            .unwrap_err()
            .detail
            .contains("axis-aligned")
    );
    invalid = original;
    invalid.atoms[0].header_extra = vec!["a".repeat(40), "b".repeat(40)];
    assert!(
        write_dlpoly_config(&invalid)
            .unwrap_err()
            .detail
            .contains("atom header exceeds 72-byte")
    );
}
