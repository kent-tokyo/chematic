use super::*;

// -----------------------------------------------------------------
// Input: round trip
// -----------------------------------------------------------------

const WATER_OPT_FREQ_INP: &str = r#"# chematic ORCA input fixture -- hand-authored per the public ORCA
# manual's documented syntax (see module docs for exact source chapters).
! B3LYP def2-SVP Opt Freq TightSCF

%pal nprocs 4 end

%scf
  MaxIter 200
  convergence tight
end

* xyz 0 1
O   0.000000   0.000000   0.000000
H   0.000000   0.757200   0.586200
H   0.000000  -0.757200   0.586200
*
"#;

#[test]
fn parse_input_keywords_concatenated_across_lines() {
    let input = "! B3LYP def2-SVP\n! Opt Freq\n\n* xyz 0 1\nO 0 0 0\n*\n";
    let parsed = parse_orca_input(input).unwrap();
    assert_eq!(parsed.keywords, vec!["B3LYP", "def2-SVP", "Opt", "Freq"]);
}

#[test]
fn parse_input_basic_fixture() {
    let parsed = parse_orca_input(WATER_OPT_FREQ_INP).unwrap();
    assert_eq!(parsed.comments.len(), 2);
    assert_eq!(
        parsed.keywords,
        vec!["B3LYP", "def2-SVP", "Opt", "Freq", "TightSCF"]
    );
    assert_eq!(parsed.blocks.len(), 2);
    assert_eq!(parsed.blocks[0].name, "pal");
    assert_eq!(parsed.blocks[0].raw, "nprocs 4");
    assert_eq!(parsed.blocks[1].name, "scf");
    assert!(parsed.blocks[1].raw.contains("MaxIter 200"));
    assert!(parsed.blocks[1].raw.contains("convergence tight"));

    let coords = parsed.coords.expect("coordinate block");
    match &coords {
        OrcaCoords::Xyz {
            charge,
            multiplicity,
            atoms,
        } => {
            assert_eq!(*charge, 0);
            assert_eq!(*multiplicity, 1);
            assert_eq!(atoms.len(), 3);
            assert_eq!(atoms[0].element, Element::from_symbol("O").unwrap());
            assert_eq!(atoms[1].element, Element::from_symbol("H").unwrap());
        }
        other => panic!("expected Xyz coords, got {other:?}"),
    }

    let (mol, xyz, charge, mult) = coords.to_molecule().expect("xyz -> molecule");
    assert_eq!(mol.atom_count(), 3);
    assert_eq!(xyz.len(), 3);
    assert_eq!(charge, 0);
    assert_eq!(mult, 1);
}

#[test]
fn input_round_trip_semantic_equality() {
    let parsed1 = parse_orca_input(WATER_OPT_FREQ_INP).unwrap();
    let written = write_orca_input(&parsed1);
    let parsed2 = parse_orca_input(&written).unwrap();
    assert_eq!(parsed1, parsed2);
}

#[test]
fn input_round_trip_unknown_block_survives_verbatim() {
    // %output is deliberately "unknown" to this parser (it doesn't
    // understand any block's internal syntax) -- it must be preserved
    // as opaque raw text, not interpreted or dropped, including a
    // value that happens to look like it could be misread (a quoted
    // string containing "end" as a substring of a longer word).
    let input = "! Opt\n\n%output\n  Print[ P_Mulliken ] 1\n  Print[ P_Loewdin ] 1\n  Appendix \"trend.txt\"\nend\n\n* xyz 0 1\nC 0 0 0\n*\n";
    let parsed1 = parse_orca_input(input).unwrap();
    assert_eq!(parsed1.blocks.len(), 1);
    assert_eq!(parsed1.blocks[0].name, "output");
    assert!(parsed1.blocks[0].raw.contains("P_Mulliken"));
    assert!(parsed1.blocks[0].raw.contains("P_Loewdin"));
    assert!(parsed1.blocks[0].raw.contains("Appendix \"trend.txt\""));

    let written = write_orca_input(&parsed1);
    let parsed2 = parse_orca_input(&written).unwrap();
    assert_eq!(parsed1.blocks, parsed2.blocks);
}

#[test]
fn input_block_with_nested_end_round_trips() {
    // Real ORCA `%geom` blocks can contain a nested `Constraints ...
    // end` sub-block that closes with its own `end` *before* the
    // block's real closing `end`. This parser finds the block's true
    // terminator as the *last* `end`-terminated line before the next
    // top-level construct, not the first -- so the inner `end` is
    // correctly captured as ordinary raw content instead of being
    // mistaken for the block's terminator (see module docs).
    let input = "! Opt\n\n%geom\n Constraints\n  {C 0 C}\n end\nend\n\n* xyz 0 1\nC 0 0 0\n*\n";
    let parsed1 = parse_orca_input(input).unwrap();
    assert_eq!(parsed1.blocks.len(), 1);
    assert_eq!(parsed1.blocks[0].name, "geom");
    assert!(parsed1.blocks[0].has_end);
    assert!(parsed1.blocks[0].raw.contains("Constraints"));
    assert!(parsed1.blocks[0].raw.contains("{C 0 C}"));
    // The inner `end` line survives verbatim as raw content.
    assert!(parsed1.blocks[0].raw.trim_end().ends_with("end"));

    // The coordinate block after it must still parse normally -- the
    // outer `end` was correctly consumed as the block's terminator,
    // not left as a stray top-level line.
    let coords = parsed1.coords.as_ref().expect("coordinate block");
    assert!(matches!(coords, OrcaCoords::Xyz { .. }));

    let written = write_orca_input(&parsed1);
    let parsed2 = parse_orca_input(&written).unwrap();
    assert_eq!(parsed1, parsed2);
}

#[test]
fn input_no_end_single_line_directive_round_trips() {
    // `%maxcore <n>` is a single-line ORCA "global variable" directive
    // with no `end` at all (unlike `%pal ... end`) -- confirmed via
    // multiple independent secondary sources, see module docs.
    let input = "! Opt\n\n%maxcore 3000\n%pal nprocs 4 end\n\n* xyz 0 1\nHe 0 0 0\n*\n";
    let parsed1 = parse_orca_input(input).unwrap();
    assert_eq!(parsed1.blocks.len(), 2);
    assert_eq!(parsed1.blocks[0].name, "maxcore");
    assert!(!parsed1.blocks[0].has_end);
    assert_eq!(parsed1.blocks[0].raw, "3000");
    assert_eq!(parsed1.blocks[1].name, "pal");
    assert!(parsed1.blocks[1].has_end);

    let written = write_orca_input(&parsed1);
    assert!(written.contains("%maxcore 3000\n"));
    assert!(!written.contains("%maxcore 3000 end"));
    let parsed2 = parse_orca_input(&written).unwrap();
    assert_eq!(parsed1, parsed2);
}

#[test]
fn input_no_end_directive_followed_by_multiline_block() {
    // Sharper version of the above: the `has_end: false` %maxcore
    // block's search window must stop at the following %scf block's
    // own line, never reach past it to grab %scf's `end` a few lines
    // later. If windowing were broken, %maxcore would swallow %scf's
    // content and %scf would be left with no closing `end` at all.
    let input = "! Opt\n\n%maxcore 3000\n%scf\n  MaxIter 200\n  convergence tight\nend\n\n* xyz 0 1\nHe 0 0 0\n*\n";
    let parsed1 = parse_orca_input(input).unwrap();
    assert_eq!(parsed1.blocks.len(), 2);
    assert_eq!(parsed1.blocks[0].name, "maxcore");
    assert!(!parsed1.blocks[0].has_end);
    assert_eq!(parsed1.blocks[0].raw, "3000");
    assert_eq!(parsed1.blocks[1].name, "scf");
    assert!(parsed1.blocks[1].has_end);
    assert!(parsed1.blocks[1].raw.contains("MaxIter 200"));
    assert!(parsed1.blocks[1].raw.contains("convergence tight"));

    let written = write_orca_input(&parsed1);
    let parsed2 = parse_orca_input(&written).unwrap();
    assert_eq!(parsed1, parsed2);
}

#[test]
fn input_round_trip_frozen_coordinate_marker() {
    let input =
        "! Opt\n\n* xyz 0 1\nC 0.000000$ 0.000000$ 0.000000$\nO 1.200000 0.000000 0.000000\n*\n";
    let parsed1 = parse_orca_input(input).unwrap();
    let OrcaCoords::Xyz { atoms, .. } = parsed1.coords.as_ref().unwrap() else {
        panic!("expected xyz coords");
    };
    assert_eq!(atoms[0].frozen, [true, true, true]);
    assert_eq!(atoms[1].frozen, [false, false, false]);

    let written = write_orca_input(&parsed1);
    assert!(written.contains("0.000000$"));
    let parsed2 = parse_orca_input(&written).unwrap();
    assert_eq!(parsed1, parsed2);
}

#[test]
fn input_single_line_block_round_trips() {
    let input = "! sp\n\n%pal nprocs 8 end\n\n* xyz 0 1\nHe 0 0 0\n*\n";
    let parsed1 = parse_orca_input(input).unwrap();
    assert_eq!(parsed1.blocks[0].raw, "nprocs 8");
    let written = write_orca_input(&parsed1);
    let parsed2 = parse_orca_input(&written).unwrap();
    assert_eq!(parsed1, parsed2);
}

#[test]
fn input_xyzfile_variant() {
    let input = "! sp\n\n* xyzfile 0 1 geom.xyz\n";
    let parsed = parse_orca_input(input).unwrap();
    match parsed.coords.unwrap() {
        OrcaCoords::XyzFile {
            charge,
            multiplicity,
            filename,
        } => {
            assert_eq!(charge, 0);
            assert_eq!(multiplicity, 1);
            assert_eq!(filename, "geom.xyz");
        }
        other => panic!("expected XyzFile, got {other:?}"),
    }
}

#[test]
fn bounded_input_parser_rejects_input_and_line_limits() {
    assert!(matches!(
        parse_orca_input_with_limits(
            WATER_OPT_FREQ_INP,
            &OrcaInputParseLimits {
                max_input_bytes: 8,
                ..Default::default()
            }
        ),
        Err(OrcaInputError::ResourceLimit {
            resource: "input bytes",
            ..
        })
    ));

    let long_line = format!("{}\n", "x".repeat(32));
    assert!(matches!(
        parse_orca_input_with_limits(
            &long_line,
            &OrcaInputParseLimits {
                max_line_bytes: 16,
                ..Default::default()
            }
        ),
        Err(OrcaInputError::ResourceLimit {
            resource: "line bytes",
            ..
        })
    ));
}

#[test]
fn bounded_input_parser_rejects_coordinate_atom_limit() {
    assert!(matches!(
        parse_orca_input_with_limits(
            WATER_OPT_FREQ_INP,
            &OrcaInputParseLimits {
                max_atoms: 2,
                ..Default::default()
            }
        ),
        Err(OrcaInputError::ResourceLimit {
            resource: "atoms",
            ..
        })
    ));
}

#[test]
fn input_int_variant_preserves_raw_zmatrix() {
    let input = "! sp\n\n* int 0 1\nO 0 0 0 0.0 0.0 0.0\nH 1 0 0 0.96 0.0 0.0\n*\n";
    let parsed = parse_orca_input(input).unwrap();
    match parsed.coords.unwrap() {
        OrcaCoords::Internal {
            charge,
            multiplicity,
            raw,
        } => {
            assert_eq!(charge, 0);
            assert_eq!(multiplicity, 1);
            assert!(raw.contains("H 1 0 0"));
        }
        other => panic!("expected Internal, got {other:?}"),
    }
}

// -----------------------------------------------------------------
// Input: malformed / doesn't panic
// -----------------------------------------------------------------

#[test]
fn input_unterminated_block_is_typed_error_not_panic() {
    let input = "! Opt\n\n%scf\n  convergence tight\n";
    let err = parse_orca_input(input).unwrap_err();
    assert_eq!(
        err,
        OrcaInputError::UnterminatedBlock {
            name: "scf".to_string()
        }
    );
}

#[test]
fn input_unterminated_coord_block_is_typed_error() {
    let input = "! Opt\n\n* xyz 0 1\nC 0 0 0\n";
    let err = parse_orca_input(input).unwrap_err();
    assert_eq!(err, OrcaInputError::UnterminatedCoordBlock);
}

#[test]
fn input_unknown_element_is_typed_error() {
    let input = "! Opt\n\n* xyz 0 1\nZq 0 0 0\n*\n";
    let err = parse_orca_input(input).unwrap_err();
    assert!(matches!(err, OrcaInputError::UnknownElement { .. }));
}

#[test]
fn input_nan_infinity_coordinate_rejected() {
    for bad in ["NaN", "inf", "-inf", "infinity"] {
        let input = format!("! Opt\n\n* xyz 0 1\nC {bad} 0.0 0.0\n*\n");
        let err = parse_orca_input(&input).unwrap_err();
        assert!(
            matches!(
                err,
                OrcaInputError::NonFiniteCoordinate { .. }
                    | OrcaInputError::InvalidCoordinate { .. }
            ),
            "unexpected error for {bad}: {err:?}"
        );
    }
}

#[test]
fn input_garbage_top_level_line_is_typed_error_not_panic() {
    let input = "this is not valid orca input at all\n";
    let err = parse_orca_input(input).unwrap_err();
    assert!(matches!(err, OrcaInputError::UnexpectedLine { .. }));
}

#[test]
fn input_empty_string_parses_to_empty_struct() {
    let parsed = parse_orca_input("").unwrap();
    assert_eq!(parsed, OrcaInput::default());
}

#[test]
fn input_malformed_charge_multiplicity_is_typed_error() {
    let input = "* xyz notanumber 1\nC 0 0 0\n*\n";
    let err = parse_orca_input(input).unwrap_err();
    assert_eq!(err, OrcaInputError::MalformedCoordHeader { line: 1 });
}

// -----------------------------------------------------------------
// Output: fixtures
// -----------------------------------------------------------------

/// Synthetic-but-realistic ORCA output for a water B3LYP/def2-SVP
/// `Opt Freq` job that terminates normally with a converged geometry.
/// Modeled on the following *confirmed* real ORCA output conventions
/// (see module docs for the exact source and confidence level of
/// each): the `Total Charge`/`Multiplicity` echo, the
/// `GEOMETRY OPTIMIZATION CYCLE N` step header, the
/// `CARTESIAN COORDINATES (ANGSTROEM)` block (verified against a real
/// reference ORCA output file, not just manual prose), the
/// `...HURRAY.../THE OPTIMIZATION HAS CONVERGED` banner, the
/// `FINAL SINGLE POINT ENERGY` line, the `VIBRATIONAL FREQUENCIES`
/// block, and the `****ORCA TERMINATED NORMALLY****` banner. The
/// energy/frequency/geometry *numeric values* are illustrative
/// (plausible water B3LYP/def2-SVP-scale numbers), not taken from an
/// actual ORCA run (no ORCA installation was available to generate
/// one) -- only the surrounding textual conventions are meant to be
/// realistic.
const WATER_OPT_FREQ_OUT: &str = r#"
                                 * O   R   C   A *
                        -- An Ab Initio, DFT and Semiempirical program --

 Total Charge           Charge          ....    0
 Multiplicity           Mult            ....    1

         *************************************************************
         *                GEOMETRY OPTIMIZATION CYCLE   1            *
         *************************************************************

---------------------------------
CARTESIAN COORDINATES (ANGSTROEM)
---------------------------------
  O     -0.000027   -0.086796   -1.479708
  H      0.000079    0.826472   -1.821557
  H     -0.000140    0.055788   -0.514893

-------------------------   --------------------
FINAL SINGLE POINT ENERGY       -76.320145981234
-------------------------   --------------------

         *************************************************************
         *                GEOMETRY OPTIMIZATION CYCLE   2            *
         *************************************************************

---------------------------------
CARTESIAN COORDINATES (ANGSTROEM)
---------------------------------
  O     -0.000015   -0.065012   -1.469881
  H      0.000041    0.798224   -1.798113
  H     -0.000090    0.071366   -0.520994

-------------------------   --------------------
FINAL SINGLE POINT ENERGY       -76.325509112233
-------------------------   --------------------

                    ***********************HURRAY********************
                    ***        THE OPTIMIZATION HAS CONVERGED     ***
                    *************************************************

---------------------------------
CARTESIAN COORDINATES (ANGSTROEM)
---------------------------------
  O     -0.000012   -0.063981   -1.468552
  H      0.000038    0.795871   -1.795220
  H     -0.000085    0.072998   -0.522004

-------------------------   --------------------
FINAL SINGLE POINT ENERGY       -76.325611004521
-------------------------   --------------------

-----------------------
VIBRATIONAL FREQUENCIES
-----------------------

   0:         0.00 cm**-1
   1:         0.00 cm**-1
   2:         0.00 cm**-1
   3:         0.00 cm**-1
   4:         0.00 cm**-1
   5:         0.00 cm**-1
   6:      1614.32 cm**-1
   7:      3672.88 cm**-1
   8:      3785.41 cm**-1

                    ****ORCA TERMINATED NORMALLY****
TOTAL RUN TIME: 0 days 0 hours 0 minutes 4 seconds 118 msec
"#;

#[test]
fn output_charge_and_multiplicity() {
    let out = parse_orca_output(WATER_OPT_FREQ_OUT).unwrap();
    assert_eq!(out.charge, Some(0));
    assert_eq!(out.multiplicity, Some(1));
}

#[test]
fn output_final_energy() {
    let out = parse_orca_output(WATER_OPT_FREQ_OUT).unwrap();
    let e = out.final_energy_hartree.expect("energy present");
    assert!((e - (-76.325611004521)).abs() < 1e-9);
}

#[test]
fn output_final_geometry_and_trajectory() {
    let out = parse_orca_output(WATER_OPT_FREQ_OUT).unwrap();
    // 3 CARTESIAN COORDINATES (ANGSTROEM) blocks in the fixture.
    assert_eq!(out.trajectory.len(), 3);
    for frame in &out.trajectory {
        assert_eq!(frame.mol.atom_count(), 3);
        assert_eq!(frame.coords.len(), 3);
    }
    let last = out.final_geometry().expect("final geometry");
    let (x, y, z) = last.coords[0];
    assert!((x - (-0.000012)).abs() < 1e-6);
    assert!((y - (-0.063981)).abs() < 1e-6);
    assert!((z - (-1.468552)).abs() < 1e-6);
}

#[test]
fn output_frequencies() {
    let out = parse_orca_output(WATER_OPT_FREQ_OUT).unwrap();
    assert_eq!(out.frequencies_cm1.len(), 9);
    assert_eq!(out.frequencies_cm1[0], 0.0);
    assert!((out.frequencies_cm1[6] - 1614.32).abs() < 1e-6);
    assert!((out.frequencies_cm1[8] - 3785.41).abs() < 1e-6);
}

#[test]
fn output_termination_and_convergence() {
    let out = parse_orca_output(WATER_OPT_FREQ_OUT).unwrap();
    assert_eq!(out.termination, OrcaTermination::Normal);
    assert_eq!(out.optimization_convergence, OrcaOptConvergence::Converged);
}

#[test]
fn output_truncated_mid_optimization_is_incomplete_not_error() {
    // Same job, but the log ends abruptly mid-SCF-cycle -- as if the
    // job crashed or was killed. No HURRAY, no termination banner.
    let cut = WATER_OPT_FREQ_OUT
        .split("GEOMETRY OPTIMIZATION CYCLE   2")
        .next()
        .unwrap();
    let out = parse_orca_output(cut).unwrap();
    assert_eq!(out.termination, OrcaTermination::Incomplete);
    assert_eq!(out.optimization_convergence, OrcaOptConvergence::Unknown);
    // The one geometry block before the cut is still recovered.
    assert_eq!(out.trajectory.len(), 1);
}

#[test]
fn output_error_termination_detected() {
    let input = "\n Total Charge           Charge          ....    0\n Multiplicity           Mult            ....    1\n\nORCA finished by error termination in SCF\nCalling Command: orca_scf job.scf.tmp job\n";
    let out = parse_orca_output(input).unwrap();
    match out.termination {
        OrcaTermination::Error(msg) => {
            assert!(msg.contains("ORCA finished by error termination"))
        }
        other => panic!("expected Error termination, got {other:?}"),
    }
}

#[test]
fn output_optimization_not_converged() {
    let input = "\n         *************************************************************\n         *                GEOMETRY OPTIMIZATION CYCLE  50            *\n         *************************************************************\n\nWarning\n   The optimization did not converge but reached the maximum number of\n   optimization cycles. Please check your results very carefully.\n\n****ORCA TERMINATED NORMALLY****\n";
    let out = parse_orca_output(input).unwrap();
    assert_eq!(
        out.optimization_convergence,
        OrcaOptConvergence::NotConverged
    );
    assert_eq!(out.termination, OrcaTermination::Normal);
}

#[test]
fn output_no_optimization_job_is_not_requested() {
    let input = "\n Total Charge           Charge          ....    0\n Multiplicity           Mult            ....    1\n\n-------------------------   --------------------\nFINAL SINGLE POINT ENERGY       -76.025678123456\n-------------------------   --------------------\n\n****ORCA TERMINATED NORMALLY****\n";
    let out = parse_orca_output(input).unwrap();
    assert_eq!(
        out.optimization_convergence,
        OrcaOptConvergence::NotRequested
    );
    assert!((out.final_energy_hartree.unwrap() - (-76.025678123456)).abs() < 1e-9);
}

#[test]
fn output_empty_input_is_incomplete_not_panic() {
    let out = parse_orca_output("").unwrap();
    assert_eq!(out.termination, OrcaTermination::Incomplete);
    assert_eq!(
        out.optimization_convergence,
        OrcaOptConvergence::NotRequested
    );
    assert!(out.trajectory.is_empty());
    assert!(out.frequencies_cm1.is_empty());
}

#[test]
fn output_garbage_input_does_not_panic() {
    let inputs = [
        "not an orca file at all",
        "CARTESIAN COORDINATES (ANGSTROEM)\n",
        "VIBRATIONAL FREQUENCIES\n-----\n",
        "* garbage * garbage *\n\x00\x01binary junk",
        "FINAL SINGLE POINT ENERGY\n",
    ];
    for input in inputs {
        let _ = parse_orca_output(input);
    }
}

#[test]
fn output_unrecognized_fragment_is_incomplete_result() {
    let out = parse_orca_output("partial launcher output without an ORCA marker").unwrap();
    assert_eq!(out.termination, OrcaTermination::Incomplete);
    assert_eq!(out.final_energy_hartree, None);
    assert!(out.trajectory.is_empty());
}

/// [`OrcaOutput`] doesn't implement `Debug` (it holds a [`Molecule`],
/// which doesn't either -- same convention as `GaussianLogResult` in
/// `gaussian.rs`), so `Result::unwrap_err` can't be used directly on
/// [`parse_orca_output`]'s return value. This is the test-only stand-in.
fn expect_output_err(r: Result<OrcaOutput, OrcaOutputError>) -> OrcaOutputError {
    match r {
        Ok(_) => panic!("expected an error, got Ok"),
        Err(e) => e,
    }
}

#[test]
fn output_nan_infinity_energy_rejected() {
    for bad in ["NaN", "inf", "-inf"] {
        let input = format!("FINAL SINGLE POINT ENERGY       {bad}\n");
        let err = expect_output_err(parse_orca_output(&input));
        assert!(matches!(err, OrcaOutputError::NonFiniteEnergy(_)));
    }
}

#[test]
fn output_nan_infinity_coordinate_rejected() {
    let input = "CARTESIAN COORDINATES (ANGSTROEM)\n---\n  O     NaN   0.0   0.0\n";
    let err = expect_output_err(parse_orca_output(input));
    assert!(matches!(err, OrcaOutputError::NonFiniteCoordinate { .. }));
}

#[test]
fn output_nan_infinity_frequency_rejected() {
    let input = "VIBRATIONAL FREQUENCIES\n-----\n\n   0:         NaN cm**-1\n";
    let err = expect_output_err(parse_orca_output(input));
    assert!(matches!(err, OrcaOutputError::NonFiniteFrequency(_)));
}

#[test]
fn output_size_limit_enforced() {
    let huge = "x\n".repeat(MAX_OUTPUT_LINES + 10);
    let err = expect_output_err(parse_orca_output(&huge));
    assert!(matches!(err, OrcaOutputError::TooManyLines { .. }));
}

#[test]
fn output_explicit_limits_reject_line_and_frame_growth() {
    let err = expect_output_err(parse_orca_output_with_limits(
        WATER_OPT_FREQ_OUT,
        &OrcaOutputParseLimits {
            max_line_bytes: 8,
            ..Default::default()
        },
    ));
    assert!(matches!(err, OrcaOutputError::LineTooLong { .. }));

    let err = expect_output_err(parse_orca_output_with_limits(
        WATER_OPT_FREQ_OUT,
        &OrcaOutputParseLimits {
            max_geometry_frames: 2,
            ..Default::default()
        },
    ));
    assert!(matches!(err, OrcaOutputError::TooManyGeometryFrames { .. }));
}

#[test]
fn output_explicit_limits_reject_geometry_and_frequency_growth() {
    let err = expect_output_err(parse_orca_output_with_limits(
        WATER_OPT_FREQ_OUT,
        &OrcaOutputParseLimits {
            max_geometry_atoms: 2,
            ..Default::default()
        },
    ));
    assert!(matches!(err, OrcaOutputError::TooManyGeometryAtoms { .. }));

    let err = expect_output_err(parse_orca_output_with_limits(
        WATER_OPT_FREQ_OUT,
        &OrcaOutputParseLimits {
            max_frequencies: 8,
            ..Default::default()
        },
    ));
    assert!(matches!(err, OrcaOutputError::TooManyFrequencies { .. }));
}
