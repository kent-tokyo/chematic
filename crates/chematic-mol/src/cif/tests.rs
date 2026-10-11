use super::*;
use crate::cif_symmetry;
use chematic_crystal::CrystalError;

const NACL_P1_CIF: &str = "data_NaCl\n\
            _cell_length_a 5.6402\n_cell_length_b 5.6402\n_cell_length_c 5.6402\n\
            _cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n\
            loop_\n_atom_site_label\n_atom_site_type_symbol\n\
            _atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n\
            Na1 Na 0.0 0.0 0.0\nCl1 Cl 0.5 0.5 0.5\n";

#[test]
fn parses_site_count_and_lattice() {
    let r = parse_cif_periodic_structure(NACL_P1_CIF).unwrap();
    assert_eq!(r.structure.site_count(), 2);
    assert!((r.structure.lattice().lengths()[0] - 5.6402).abs() < 1e-4);
    assert_eq!(r.symmetry, CifSymmetryStatus::P1);
}

#[test]
fn no_symmetry_tags_at_all_classifies_as_p1() {
    // NACL_P1_CIF has zero _symmetry_*/_space_group_* tags -- matches
    // the existing parser's "effectively P1" scope for such files.
    let r = parse_cif_periodic_structure(NACL_P1_CIF).unwrap();
    assert_eq!(r.symmetry, CifSymmetryStatus::P1);
}

/// Shared C2/c (space group No. 15) fixture: a realistic 8-operation
/// symop list with one asymmetric-unit atom sitting on a special
/// position. Used both to prove the default (expand=true) path
/// genuinely expands+dedups, and (via `expand_explicit_symmetry:
/// false`) that opting out still reports the pre-expansion
/// asymmetric-unit-only behavior.
const C2C_CIF: &str = "data_synthetic_c2c\n\
            _cell_length_a 10.0\n_cell_length_b 8.0\n_cell_length_c 12.0\n\
            _cell_angle_alpha 90\n_cell_angle_beta 105.0\n_cell_angle_gamma 90\n\
            _symmetry_space_group_name_H-M 'C 2/c'\n\
            loop_\n_space_group_symop_operation_xyz\n\
            'x, y, z'\n'-x, y, -z+1/2'\n'-x, -y, -z'\n'x, -y, z+1/2'\n\
            'x+1/2, y+1/2, z'\n'-x+1/2, y+1/2, -z+1/2'\n'-x+1/2, -y+1/2, -z'\n\
            'x+1/2, -y+1/2, z+1/2'\n\
            loop_\n_atom_site_label\n_atom_site_type_symbol\n\
            _atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n\
            Ti1 Ti 0.0 0.25 0.25\n";

#[test]
fn single_identity_operation_cif_expands_trivially_by_default() {
    // A genuinely P1-equivalent CIF that explicitly lists its one
    // (identity) operation -- IUCr's own convention for how a P1
    // file states its symmetry explicitly. Under the default
    // (expand_explicit_symmetry=true), this is not P1-classified by
    // absence of an operation list; it is honestly reported as
    // ExpandedExplicitOperations, expanding trivially to itself
    // (never an error -- see CifSymmetryError::MissingIdentityOperation's
    // docs for why a *present* identity, even alone, is not an error).
    let cif = "data_x\n\
                _cell_length_a 5.0\n_cell_length_b 5.0\n_cell_length_c 5.0\n\
                _cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n\
                _symmetry_space_group_name_H-M 'P 1'\n\
                loop_\n_space_group_symop_operation_xyz\n'x, y, z'\n\
                loop_\n_atom_site_label\n_atom_site_type_symbol\n\
                _atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n\
                C1 C 0.0 0.0 0.0\n";
    let r = parse_cif_periodic_structure(cif).unwrap();
    assert_eq!(r.structure.site_count(), 1);
    assert_eq!(
        r.symmetry,
        CifSymmetryStatus::ExpandedExplicitOperations {
            space_group_name: Some("P 1".to_string()),
            operation_count: 1,
            asymmetric_site_count: 1,
            expanded_site_count: 1,
        }
    );
    assert!(r.symmetry.is_complete_cell());
}

#[test]
fn single_identity_operation_cif_reports_unexpanded_when_expansion_is_opted_out() {
    let cif = "data_x\n\
                _cell_length_a 5.0\n_cell_length_b 5.0\n_cell_length_c 5.0\n\
                _cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n\
                _symmetry_space_group_name_H-M 'P 1'\n\
                loop_\n_space_group_symop_operation_xyz\n'x, y, z'\n\
                loop_\n_atom_site_label\n_atom_site_type_symbol\n\
                _atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n\
                C1 C 0.0 0.0 0.0\n";
    let r = parse_cif_periodic_structure_with_options(
        cif,
        CifPeriodicParseOptions {
            expand_explicit_symmetry: false,
        },
    )
    .unwrap();
    assert_eq!(r.structure.site_count(), 1);
    assert_eq!(
        r.symmetry,
        CifSymmetryStatus::UnexpandedSymmetry {
            space_group_name: Some("P 1".to_string()),
            operation_count: 1,
        }
    );
    assert!(!r.symmetry.is_complete_cell());
}

#[test]
fn declared_non_p1_space_group_expands_by_default_and_dedups_a_special_position() {
    let r = parse_cif_periodic_structure(C2C_CIF).unwrap();
    // Ti1 sits on a special position of C2/c (Wyckoff 4e): the
    // 8-operation general-position list collapses pairwise, so the
    // expanded cell has 4 distinct sites, not 8.
    assert_eq!(r.structure.site_count(), 4);
    match r.symmetry {
        CifSymmetryStatus::ExpandedExplicitOperations {
            space_group_name,
            operation_count,
            asymmetric_site_count,
            expanded_site_count,
        } => {
            assert_eq!(space_group_name.as_deref(), Some("C 2/c"));
            assert_eq!(operation_count, 8);
            assert_eq!(asymmetric_site_count, 1);
            assert_eq!(expanded_site_count, 4);
        }
        other => panic!("expected ExpandedExplicitOperations, got {other:?}"),
    }
}

#[test]
fn declared_non_p1_space_group_is_not_silently_treated_as_p1_when_not_expanded() {
    let r = parse_cif_periodic_structure_with_options(
        C2C_CIF,
        CifPeriodicParseOptions {
            expand_explicit_symmetry: false,
        },
    )
    .unwrap();
    assert_eq!(r.structure.site_count(), 1);
    match r.symmetry {
        CifSymmetryStatus::UnexpandedSymmetry {
            space_group_name,
            operation_count,
        } => {
            assert_eq!(space_group_name.as_deref(), Some("C 2/c"));
            assert_eq!(operation_count, 8);
        }
        other => panic!("expected UnexpandedSymmetry, got {other:?}"),
    }
}

#[test]
fn it_number_alone_flags_non_p1_even_without_a_name_or_symop_loop() {
    let cif = "data_x\n\
                _cell_length_a 5.0\n_cell_length_b 5.0\n_cell_length_c 5.0\n\
                _cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n\
                _symmetry_Int_Tables_number 15\n\
                loop_\n_atom_site_label\n_atom_site_type_symbol\n\
                _atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n\
                Ti1 Ti 0.0 0.25 0.25\n";
    let r = parse_cif_periodic_structure(cif).unwrap();
    match r.symmetry {
        CifSymmetryStatus::UnexpandedSymmetry {
            space_group_name,
            operation_count,
        } => {
            assert_eq!(space_group_name, None);
            assert_eq!(operation_count, 0);
        }
        other => panic!("expected UnexpandedSymmetry, got {other:?}"),
    }
}

#[test]
fn disorder_rows_at_the_same_site_merge_into_multi_species_site() {
    // Realistic disorder loop shape: extra columns
    // (U_iso_or_equiv/adp_type/symmetry_multiplicity/calc_flag)
    // interleaved before occupancy, as commonly exported by COD/ICSD
    // -- exercises column-offset resolution, not just a 5-column toy.
    let cif = "data_synthetic_disordered_alloy\n\
                _cell_length_a 4.0\n_cell_length_b 4.0\n_cell_length_c 4.0\n\
                _cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n\
                loop_\n\
                _atom_site_label\n\
                _atom_site_type_symbol\n\
                _atom_site_fract_x\n\
                _atom_site_fract_y\n\
                _atom_site_fract_z\n\
                _atom_site_U_iso_or_equiv\n\
                _atom_site_adp_type\n\
                _atom_site_occupancy\n\
                _atom_site_symmetry_multiplicity\n\
                _atom_site_calc_flag\n\
                Fe1 Fe 0.0 0.0 0.0 0.012 Uiso 0.6 1 d\n\
                Ni1 Ni 0.0 0.0 0.0 0.012 Uiso 0.4 1 d\n\
                O1  O  0.5 0.5 0.5 0.020 Uiso 1.0 1 d\n";
    let r = parse_cif_periodic_structure(cif).unwrap();
    assert_eq!(
        r.structure.site_count(),
        2,
        "Fe1/Ni1 must merge into one site"
    );
    let disordered = &r.structure.sites()[0];
    assert_eq!(disordered.species.len(), 2);
    assert_eq!(disordered.label.as_deref(), Some("Fe1"));
    assert_eq!(disordered.species[0].element, Element::FE);
    assert_eq!(disordered.species[1].element, Element::NI);
    let occs: Vec<f64> = disordered
        .species
        .iter()
        .map(|s| s.occupancy.value())
        .collect();
    assert!((occs[0] - 0.6).abs() < 1e-9);
    assert!((occs[1] - 0.4).abs() < 1e-9);
    let ordered = &r.structure.sites()[1];
    assert_eq!(ordered.species.len(), 1);
    assert_eq!(ordered.species[0].element, Element::O);
    assert!((ordered.species[0].occupancy.value() - 1.0).abs() < 1e-9);
}

#[test]
fn occupancy_sum_exceeded_is_rejected_by_crystals_own_validator() {
    let cif = "data_x\n\
                _cell_length_a 4.0\n_cell_length_b 4.0\n_cell_length_c 4.0\n\
                _cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n\
                loop_\n_atom_site_label\n_atom_site_type_symbol\n\
                _atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n\
                _atom_site_occupancy\n\
                Fe1 Fe 0.0 0.0 0.0 0.7\nNi1 Ni 0.0 0.0 0.0 0.5\n";
    let err = parse_cif_periodic_structure(cif).unwrap_err();
    // Raised by PeriodicSite::new itself (before PeriodicStructure
    // wraps per-site errors in InvalidSite{index,..}), proving this
    // adapter reuses chematic-crystal's own occupancy-sum validator
    // rather than re-checking the sum a second, different way.
    match err {
        CifPeriodicError::Crystal(CrystalError::OccupancySumExceeded { sum, .. }) => {
            assert!((sum - 1.2).abs() < 1e-9);
        }
        other => panic!("expected Crystal(OccupancySumExceeded), got {other:?}"),
    }
}

#[test]
fn occupancy_placeholder_dot_and_question_mark_default_to_full() {
    let cif = "data_x\n\
                _cell_length_a 4.0\n_cell_length_b 4.0\n_cell_length_c 4.0\n\
                _cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n\
                loop_\n_atom_site_label\n_atom_site_type_symbol\n\
                _atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n\
                _atom_site_occupancy\n\
                C1 C 0.0 0.0 0.0 .\nC2 C 0.5 0.5 0.5 ?\n";
    let r = parse_cif_periodic_structure(cif).unwrap();
    for site in r.structure.sites() {
        assert!((site.species[0].occupancy.value() - 1.0).abs() < 1e-9);
    }
}

#[test]
fn malformed_occupancy_value_is_an_explicit_error_not_a_silent_default() {
    let cif = "data_x\n\
                _cell_length_a 4.0\n_cell_length_b 4.0\n_cell_length_c 4.0\n\
                _cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n\
                loop_\n_atom_site_label\n_atom_site_type_symbol\n\
                _atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n\
                _atom_site_occupancy\n\
                C1 C 0.0 0.0 0.0 0.5x\n";
    let err = parse_cif_periodic_structure(cif).unwrap_err();
    assert_eq!(err, CifPeriodicError::InvalidOccupancy("0.5x".to_string()));
}

#[test]
fn missing_one_of_six_cell_tags_is_rejected_even_though_parse_cif_defaults_it() {
    // parse_cif would silently default the missing beta to 90.0;
    // the periodic adapter must not build a Lattice from a
    // partially-defaulted cell.
    let cif = "data_x\n\
                _cell_length_a 5.0\n_cell_length_b 5.0\n_cell_length_c 5.0\n\
                _cell_angle_alpha 90\n_cell_angle_gamma 90\n\
                loop_\n_atom_site_label\n_atom_site_type_symbol\n\
                _atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n\
                C1 C 0.0 0.0 0.0\n";
    let err = parse_cif_periodic_structure(cif).unwrap_err();
    assert_eq!(err, CifPeriodicError::Cif(CifError::MissingCellParameters));
    // Sanity: parse_cif (the existing, unaffected API) does NOT
    // reject this file -- confirms the two functions genuinely
    // have different cell-completeness requirements.
    assert!(super::parse_cif(cif).is_ok());
}

#[test]
fn cartesian_only_cif_converts_through_the_lattice() {
    let cif = "data_x\n\
                _cell_length_a 10.0\n_cell_length_b 10.0\n_cell_length_c 10.0\n\
                _cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n\
                loop_\n_atom_site_label\n_atom_site_type_symbol\n\
                _atom_site_Cartn_x\n_atom_site_Cartn_y\n_atom_site_Cartn_z\n\
                C1 C 5.0 5.0 5.0\n";
    let r = parse_cif_periodic_structure(cif).unwrap();
    let f = r.structure.sites()[0].fractional.0;
    assert!((f[0] - 0.5).abs() < 1e-9);
    assert!((f[1] - 0.5).abs() < 1e-9);
    assert!((f[2] - 0.5).abs() < 1e-9);
}

#[test]
fn write_then_parse_round_trips_as_p1_and_preserves_disorder() {
    let cif = "data_synthetic_disordered_alloy\n\
                _cell_length_a 4.0\n_cell_length_b 4.0\n_cell_length_c 4.0\n\
                _cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n\
                loop_\n_atom_site_label\n_atom_site_type_symbol\n\
                _atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n\
                _atom_site_occupancy\n\
                Fe1 Fe 0.0 0.0 0.0 0.6\nNi1 Ni 0.0 0.0 0.0 0.4\n\
                O1  O  0.5 0.5 0.5 1.0\n";
    let parsed = parse_cif_periodic_structure(cif).unwrap();
    let written = write_cif_periodic_structure(&parsed.structure);
    assert!(written.contains("_atom_site_occupancy"));
    let roundtripped = parse_cif_periodic_structure(&written).unwrap();
    assert_eq!(roundtripped.symmetry, CifSymmetryStatus::P1);
    assert_eq!(roundtripped.structure.site_count(), 2);
    let disordered = &roundtripped.structure.sites()[0];
    assert_eq!(disordered.species.len(), 2);
}

// -----------------------------------------------------------------
// Explicit symmetry expansion -- end-to-end (full CIF text) fixtures.
// Parser/expansion internals are exercised in isolation by
// `crate::cif_symmetry::tests`; these confirm the whole pipeline
// (tokenize -> resolve operations -> expand -> PeriodicStructure)
// agrees, using real/citable space-group data where practical.
// -----------------------------------------------------------------

fn cif_with_ops(ops_loop: &str, atom_line: &str) -> String {
    format!(
        "data_x\n\
                 _cell_length_a 6.0\n_cell_length_b 8.0\n_cell_length_c 10.0\n\
                 _cell_angle_alpha 90\n_cell_angle_beta 100.0\n_cell_angle_gamma 90\n\
                 {ops_loop}\n\
                 loop_\n_atom_site_label\n_atom_site_type_symbol\n\
                 _atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n\
                 {atom_line}\n"
    )
}

#[test]
fn p1_single_identity_operation_expands_to_itself() {
    let cif = cif_with_ops(
        "loop_\n_symmetry_equiv_pos_as_xyz\n'x,y,z'",
        "C1 C 0.1 0.2 0.3",
    );
    let r = parse_cif_periodic_structure(&cif).unwrap();
    assert_eq!(r.structure.site_count(), 1);
    assert!(matches!(
        r.symmetry,
        CifSymmetryStatus::ExpandedExplicitOperations {
            operation_count: 1,
            asymmetric_site_count: 1,
            expanded_site_count: 1,
            ..
        }
    ));
}

#[test]
fn p_minus_1_general_position_expands_to_two_sites_identity_first() {
    let cif = cif_with_ops(
        "loop_\n_symmetry_equiv_pos_as_xyz\n'x,y,z'\n'-x,-y,-z'",
        "C1 C 0.1 0.2 0.3",
    );
    let r = parse_cif_periodic_structure(&cif).unwrap();
    assert_eq!(r.structure.site_count(), 2);
    let sites = r.structure.sites();
    assert_eq!(sites[0].label.as_deref(), Some("C1"));
    let f0 = sites[0].fractional.0;
    assert!((f0[0] - 0.1).abs() < 1e-9 && (f0[1] - 0.2).abs() < 1e-9 && (f0[2] - 0.3).abs() < 1e-9);
    assert_eq!(sites[1].label.as_deref(), Some("C1@sym2"));
    let f1 = sites[1].fractional.0;
    assert!((f1[0] - 0.9).abs() < 1e-9 && (f1[1] - 0.8).abs() < 1e-9 && (f1[2] - 0.7).abs() < 1e-9);
}

#[test]
fn p_minus_1_special_position_at_origin_is_not_duplicated() {
    let cif = cif_with_ops(
        "loop_\n_symmetry_equiv_pos_as_xyz\n'x,y,z'\n'-x,-y,-z'",
        "C1 C 0.0 0.0 0.0",
    );
    let r = parse_cif_periodic_structure(&cif).unwrap();
    assert_eq!(r.structure.site_count(), 1);
    assert_eq!(r.structure.sites()[0].label.as_deref(), Some("C1"));
}

#[test]
fn cell_wrapping_translation_operation_wraps_into_unit_cell() {
    let cif = cif_with_ops(
        "loop_\n_symmetry_equiv_pos_as_xyz\n'x,y,z'\n'x+1/2,y,z'",
        "C1 C 0.8 0.1 0.1",
    );
    let r = parse_cif_periodic_structure(&cif).unwrap();
    assert_eq!(r.structure.site_count(), 2);
    let f1 = r.structure.sites()[1].fractional.0;
    assert!(
        (f1[0] - 0.3).abs() < 1e-9,
        "expected wrap to 0.3, got {}",
        f1[0]
    );
}

#[test]
fn non_trivial_rational_translation_and_variable_combination() {
    // "-y+x,-y,1/3+z": x' = x-y, y' = -y, z' = z+1/3.
    let cif = cif_with_ops(
        "loop_\n_symmetry_equiv_pos_as_xyz\n'x,y,z'\n'-y+x,-y,1/3+z'",
        "C1 C 0.5 0.6 0.7",
    );
    let r = parse_cif_periodic_structure(&cif).unwrap();
    assert_eq!(r.structure.site_count(), 2);
    let f1 = r.structure.sites()[1].fractional.0;
    assert!((f1[0] - (0.5 - 0.6f64).rem_euclid(1.0)).abs() < 1e-9);
    assert!((f1[1] - (-0.6f64).rem_euclid(1.0)).abs() < 1e-9);
    assert!((f1[2] - (0.7f64 + 1.0 / 3.0).rem_euclid(1.0)).abs() < 1e-9);
}

#[test]
fn modern_dotted_operation_tag_is_supported() {
    let cif = cif_with_ops(
        "loop_\n_space_group_symop.operation_xyz\n'x,y,z'\n'-x,-y,-z'",
        "C1 C 0.1 0.2 0.3",
    );
    let r = parse_cif_periodic_structure(&cif).unwrap();
    assert_eq!(r.structure.site_count(), 2);
}

#[test]
fn modern_dotted_equiv_pos_tag_is_supported() {
    let cif = cif_with_ops(
        "loop_\n_symmetry_equiv.pos_as_xyz\n'x,y,z'\n'-x,-y,-z'",
        "C1 C 0.1 0.2 0.3",
    );
    let r = parse_cif_periodic_structure(&cif).unwrap();
    assert_eq!(r.structure.site_count(), 2);
}

#[test]
fn legacy_underscore_equiv_pos_tag_is_supported() {
    let cif = cif_with_ops(
        "loop_\n_symmetry_equiv_pos_as_xyz\n'x,y,z'\n'-x,-y,-z'",
        "C1 C 0.1 0.2 0.3",
    );
    let r = parse_cif_periodic_structure(&cif).unwrap();
    assert_eq!(r.structure.site_count(), 2);
}

#[test]
fn trigonal_p3_three_operation_general_position_expands_to_three_sites() {
    // P3 (No. 143), standard general positions -- a real 3-operation
    // (3-fold rotation) case, distinct from the 2- and 4-operation
    // fixtures above.
    let cif = cif_with_ops(
        "loop_\n_symmetry_equiv_pos_as_xyz\n'x,y,z'\n'-y,x-y,z'\n'y-x,-x,z'",
        "C1 C 0.2 0.35 0.4",
    );
    let r = parse_cif_periodic_structure(&cif).unwrap();
    assert_eq!(r.structure.site_count(), 3);
    let labels: Vec<Option<String>> = r
        .structure
        .sites()
        .iter()
        .map(|s| s.label.clone())
        .collect();
    assert_eq!(
        labels,
        vec![
            Some("C1".to_string()),
            Some("C1@sym2".to_string()),
            Some("C1@sym3".to_string()),
        ]
    );
}

#[test]
fn fe_ni_disorder_survives_symmetry_expansion_intact() {
    // cif_with_ops's atom loop header has no occupancy column, so
    // this CIF is built directly instead of via that helper.
    let cif = "data_x\n\
                _cell_length_a 6.0\n_cell_length_b 8.0\n_cell_length_c 10.0\n\
                _cell_angle_alpha 90\n_cell_angle_beta 100.0\n_cell_angle_gamma 90\n\
                loop_\n_symmetry_equiv_pos_as_xyz\n'x,y,z'\n'-x,-y,-z'\n\
                loop_\n_atom_site_label\n_atom_site_type_symbol\n\
                _atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n\
                _atom_site_occupancy\n\
                Fe1 Fe 0.1 0.2 0.3 0.6\nNi1 Ni 0.1 0.2 0.3 0.4\n";
    let r = parse_cif_periodic_structure(cif).unwrap();
    assert_eq!(
        r.structure.site_count(),
        2,
        "one disordered site, two images"
    );
    for site in r.structure.sites() {
        assert_eq!(site.species.len(), 2);
        let fe = site
            .species
            .iter()
            .find(|s| s.element.symbol() == "Fe")
            .unwrap();
        let ni = site
            .species
            .iter()
            .find(|s| s.element.symbol() == "Ni")
            .unwrap();
        assert!((fe.occupancy.value() - 0.6).abs() < 1e-9);
        assert!((ni.occupancy.value() - 0.4).abs() < 1e-9);
    }
}

#[test]
fn duplicate_operation_is_a_typed_error() {
    let cif = cif_with_ops(
        "loop_\n_symmetry_equiv_pos_as_xyz\n'x,y,z'\n'-x,-y,-z'\n'-x,-y,-z'",
        "C1 C 0.1 0.2 0.3",
    );
    let err = parse_cif_periodic_structure(&cif).unwrap_err();
    assert!(matches!(
        err,
        CifPeriodicError::Symmetry(
            cif_symmetry::CifSymmetryError::DuplicateSymmetryOperation { .. }
        )
    ));
}

#[test]
fn missing_identity_is_a_typed_error() {
    let cif = cif_with_ops(
        "loop_\n_symmetry_equiv_pos_as_xyz\n'-x,-y,-z'",
        "C1 C 0.1 0.2 0.3",
    );
    let err = parse_cif_periodic_structure(&cif).unwrap_err();
    assert_eq!(
        err,
        CifPeriodicError::Symmetry(cif_symmetry::CifSymmetryError::MissingIdentityOperation)
    );
}

#[test]
fn declared_space_group_with_no_operation_list_at_all_stays_unexpanded() {
    // Only a name/number, no _space_group_symop*/_symmetry_equiv*
    // loop -- nothing to expand from, so this must stay
    // UnexpandedSymmetry even under the default expand=true.
    let cif = "data_x\n\
                _cell_length_a 5.0\n_cell_length_b 5.0\n_cell_length_c 5.0\n\
                _cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n\
                _symmetry_space_group_name_H-M 'P 21/c'\n\
                loop_\n_atom_site_label\n_atom_site_type_symbol\n\
                _atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n\
                C1 C 0.1 0.2 0.3\n";
    let r = parse_cif_periodic_structure(cif).unwrap();
    assert_eq!(r.structure.site_count(), 1);
    match r.symmetry {
        CifSymmetryStatus::UnexpandedSymmetry {
            space_group_name,
            operation_count,
        } => {
            assert_eq!(space_group_name.as_deref(), Some("P 21/c"));
            assert_eq!(operation_count, 0);
        }
        other => panic!("expected UnexpandedSymmetry, got {other:?}"),
    }
}

/// IUCr International Tables for Crystallography Vol. A, space
/// group P2_1/c (No. 14), standard-setting general positions --
/// public crystallographic convention, not copied from any tool's
/// implementation. Self-verifying: composes every pair of the 4
/// operations using test-local arithmetic (independent of the
/// parser under test) and confirms the result matches one of the 4
/// operations modulo integer translation, i.e. these 4 operators
/// really do form a closed group before trusting them as a golden
/// fixture.
#[test]
fn iucr_p21c_four_operation_golden_fixture() {
    type Mat = [[f64; 3]; 3];
    let identity: Mat = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let screw: Mat = [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]];
    let inversion: Mat = [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]];
    let glide: Mat = [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]];
    let ops: [(Mat, [f64; 3]); 4] = [
        (identity, [0.0, 0.0, 0.0]),
        (screw, [0.0, 0.5, 0.5]),
        (inversion, [0.0, 0.0, 0.0]),
        (glide, [0.0, 0.5, 0.5]),
    ];

    fn apply(rot: Mat, trans: [f64; 3], p: [f64; 3]) -> [f64; 3] {
        let mut out = [0.0; 3];
        for i in 0..3 {
            out[i] = trans[i] + (0..3).map(|k| rot[i][k] * p[k]).sum::<f64>();
        }
        out
    }
    fn compose(a: (Mat, [f64; 3]), b: (Mat, [f64; 3])) -> (Mat, [f64; 3]) {
        let (rot_a, t_a) = a;
        let (rot_b, t_b) = b;
        let mut rot = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                rot[i][j] = (0..3).map(|k| rot_a[i][k] * rot_b[k][j]).sum();
            }
        }
        let t = apply(rot_a, t_a, t_b);
        (rot, t)
    }
    fn mats_eq(a: Mat, b: Mat) -> bool {
        (0..3).all(|i| (0..3).all(|j| (a[i][j] - b[i][j]).abs() < 1e-9))
    }
    fn trans_eq_mod1(a: [f64; 3], b: [f64; 3]) -> bool {
        (0..3).all(|i| {
            let d = (a[i] - b[i]).rem_euclid(1.0);
            d < 1e-9 || (1.0 - d) < 1e-9
        })
    }

    // Group closure: every pairwise composition matches one of the
    // 4 declared operations, modulo integer translation.
    for &a in &ops {
        for &b in &ops {
            let (rot_c, t_c) = compose(a, b);
            let matches_some_op = ops
                .iter()
                .any(|&(rot, t)| mats_eq(rot_c, rot) && trans_eq_mod1(t_c, t));
            assert!(
                matches_some_op,
                "composition not closed within the declared set"
            );
        }
    }

    // Now confirm the actual parser+expander agrees with this
    // independently-computed reference on a general position.
    let cif = cif_with_ops(
        "loop_\n_space_group_symop_id\n_space_group_symop_operation_xyz\n\
                 1 'x,y,z'\n2 '-x,y+1/2,-z+1/2'\n3 '-x,-y,-z'\n4 'x,-y+1/2,z+1/2'",
        "C1 C 0.2 0.3 0.4",
    );
    let r = parse_cif_periodic_structure(&cif).unwrap();
    assert_eq!(r.structure.site_count(), 4);
    assert_eq!(
        r.symmetry,
        CifSymmetryStatus::ExpandedExplicitOperations {
            space_group_name: None,
            operation_count: 4,
            asymmetric_site_count: 1,
            expanded_site_count: 4,
        }
    );
    let p0 = [0.2, 0.3, 0.4];
    let expected_labels = ["C1", "C1@sym2", "C1@sym3", "C1@sym4"];
    for (i, site) in r.structure.sites().iter().enumerate() {
        let (rot, t) = ops[i];
        let expected = apply(rot, t, p0).map(|c| c.rem_euclid(1.0));
        let actual = site.fractional.0;
        for k in 0..3 {
            assert!(
                (actual[k] - expected[k]).abs() < 1e-9,
                "op {i}: expected {expected:?}, got {actual:?}"
            );
        }
        assert_eq!(site.label.as_deref(), Some(expected_labels[i]));
    }
}

// -----------------------------------------------------------------
// to_cif_checked
// -----------------------------------------------------------------

#[test]
fn to_cif_checked_succeeds_on_expanded_structure() {
    let cif = cif_with_ops(
        "loop_\n_symmetry_equiv_pos_as_xyz\n'x,y,z'\n'-x,-y,-z'",
        "C1 C 0.1 0.2 0.3",
    );
    let r = parse_cif_periodic_structure(&cif).unwrap();
    let written = r.to_cif_checked().unwrap();
    assert!(written.contains("_symmetry_space_group_name_H-M"));
    let reparsed = parse_cif_periodic_structure(&written).unwrap();
    assert_eq!(reparsed.structure.site_count(), 2);
    assert_eq!(reparsed.symmetry, CifSymmetryStatus::P1);
}

#[test]
fn to_cif_checked_refuses_unexpanded_structure() {
    let r = parse_cif_periodic_structure_with_options(
        C2C_CIF,
        CifPeriodicParseOptions {
            expand_explicit_symmetry: false,
        },
    )
    .unwrap();
    let err = r.to_cif_checked().unwrap_err();
    assert!(matches!(
        err,
        CifPeriodicError::UnexpandedSymmetryWrite {
            operation_count: 8,
            ..
        }
    ));
}

#[test]
fn to_cif_checked_succeeds_on_plain_p1() {
    let r = parse_cif_periodic_structure(NACL_P1_CIF).unwrap();
    assert!(r.to_cif_checked().is_ok());
}
