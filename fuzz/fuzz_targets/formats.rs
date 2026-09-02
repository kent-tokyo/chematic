#![no_main]

//! Dispatch one bounded input across every text parser that is cheap to invoke
//! from a standalone fuzz binary.  The selector byte keeps each execution on a
//! single parser while the corpus as a whole covers the public format surface.

use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let Some((&selector, payload)) = data.split_first() else {
        return;
    };
    let Ok(input) = std::str::from_utf8(payload) else {
        return;
    };
    let input = &input[..input.len().min(1 << 20)];

    match selector % 22 {
        0 => { let _ = chematic_chem::parse_formula(input); }
        1 => { let _ = chematic_chem::parse_condensed(input); }
        2 => { let _ = chematic_smarts::parse_smarts(input); }
        3 => { let _ = chematic_smarts::parse_cxsmarts(input); }
        4 => { let _ = chematic_rxn::parse_reaction(input); }
        5 => { let _ = chematic_rxn::parse_reaction_query(input); }
        6 => { let _ = chematic_mol::parse_cml(input); }
        7 => { let _ = chematic_mol::parse_cif(input); }
        8 => { let _ = chematic_mol::parse_cjson(input); }
        9 => { let _ = chematic_mol::parse_cube(input); }
        10 => { let _ = chematic_mol::parse_ket(input); }
        11 => { let _ = chematic_mol::parse_pdbqt(input); }
        12 => { let _ = chematic_mol::parse_pqr(input); }
        13 => { let _ = chematic_mol::parse_gjf(input); }
        14 => { let _ = chematic_mol::parse_gaussian_log(input); }
        15 => { let _ = chematic_mol::parse_mmcif(input); }
        16 => { let _ = chematic_mol::parse_mol2(input); }
        17 => { let _ = chematic_mol::parse_mrv(input); }
        18 => { let _ = chematic_mol::parse_opendx(input); }
        19 => { let _ = chematic_mol::parse_orca_input(input); }
        20 => { let _ = chematic_mol::parse_orca_output(input); }
        21 => {
            let _ = chematic_3d::parse_pdb_atoms_with_limits(
                input,
                &chematic_3d::PdbParseLimits {
                    max_input_bytes: 1 << 20,
                    max_atoms: 10_000,
                    max_line_bytes: 4096,
                    max_models: 256,
                },
            );
        }
        _ => unreachable!(),
    }
});
