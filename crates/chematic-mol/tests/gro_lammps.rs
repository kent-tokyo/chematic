use chematic_mol::gro_lammps::{GroLammpsUnits, gro_to_lammps_atomic};
use chematic_mol::{LammpsAtomStyle, parse_lammps_data, write_lammps_data};

fn row(id: u32, velocity: bool) -> String {
    let mut result = format!(
        "{:5}{:<5}{:>5}{:5}{:8.3}{:8.3}{:8.3}",
        7, "ION", "NA", id, -0.25, 0.5, 1.75
    );
    if velocity {
        result.push_str(&format!("{:8.4}{:8.4}{:8.4}", -0.125, 0.25, 0.5));
    }
    result
}
fn frame(box_line: &str, velocities: [bool; 2]) -> String {
    format!(
        "original geometry\n2\n{}\n{}\n{box_line}\n",
        row(0, velocities[0]),
        row(0, velocities[1])
    )
}

#[test]
fn orthogonal_geometry_and_velocities_use_consistent_units() {
    let text = frame("2 3 4", [true; 2]);
    for (units, factor) in [(GroLammpsUnits::Real, 0.01), (GroLammpsUnits::Metal, 10.0)] {
        let report = gro_to_lammps_atomic(&text, &[1, 2], &[22.99, 35.45], units).unwrap();
        assert_eq!(report.data.simulation_box.lo, [0.0; 3]);
        assert_eq!(report.data.simulation_box.hi, [20.0, 30.0, 40.0]);
        assert_eq!(report.data.simulation_box.tilt, None);
        assert_eq!(
            [
                report.data.atoms[0].x,
                report.data.atoms[0].y,
                report.data.atoms[0].z
            ],
            [-2.5, 5.0, 17.5]
        );
        assert_eq!(report.data.velocities[0].vx, -0.125 * factor);
        assert_eq!(report.data.velocities[0].vy, 0.25 * factor);
        assert_eq!(report.data.velocities[0].vz, 0.5 * factor);
        assert_eq!(report.data.atoms[1].atom_type, 2);
        assert_eq!(report.data.masses[1].mass, 35.45);
        assert_eq!(report.source_frame.source_text, text);
        assert_eq!(report.source_frame.atoms[0].atom_number, 0);
        assert_eq!(report.data.atoms[0].id, 1);
        assert_eq!(report.data.atoms[1].id, 2);
        let read =
            parse_lammps_data(&write_lammps_data(&report.data), LammpsAtomStyle::Atomic).unwrap();
        assert_eq!(read, report.data);
    }
}

#[test]
fn triclinic_box_retains_true_edges_and_tilt_not_enclosing_bounds() {
    let text = frame("2 3 4 0 0 0.5 0 -0.25 0.75", [false; 2]);
    let report = gro_to_lammps_atomic(&text, &[1, 1], &[12.0], GroLammpsUnits::Real).unwrap();
    assert_eq!(report.data.simulation_box.hi, [20.0, 30.0, 40.0]);
    assert_eq!(report.data.simulation_box.tilt, Some([5.0, -2.5, 7.5]));
    assert!(report.data.velocities.is_empty());
    assert!(report.data.bonds.is_empty());
    assert!(
        report
            .data
            .atoms
            .iter()
            .all(|a| a.charge.is_none() && a.molecule_id.is_none() && a.image.is_none())
    );
    let read =
        parse_lammps_data(&write_lammps_data(&report.data), LammpsAtomStyle::Atomic).unwrap();
    assert_eq!(read, report.data);
}

#[test]
fn invalid_or_unsupported_boxes_are_refused() {
    for cell in [
        "0 3 4",
        "2 -3 4",
        "2 3 0",
        "2 3 4 0.1 0 0 0 0 0",
        "2 3 4 0 0.1 0 0 0 0",
        "2 3 4 0 0 0 0.1 0 0",
        "1e308 3 4",
        "nan 3 4",
    ] {
        assert!(
            gro_to_lammps_atomic(
                &frame(cell, [false; 2]),
                &[1, 1],
                &[12.0],
                GroLammpsUnits::Metal
            )
            .is_err(),
            "{cell}"
        );
    }
}

#[test]
fn assignments_masses_and_mixed_velocities_are_refused() {
    let text = frame("2 3 4", [false; 2]);
    for types in [
        vec![],
        vec![1],
        vec![0, 1],
        vec![-1, 1],
        vec![2, 1],
        vec![i64::MAX, 1],
    ] {
        assert!(gro_to_lammps_atomic(&text, &types, &[12.0], GroLammpsUnits::Real).is_err());
    }
    for masses in [
        vec![],
        vec![0.0],
        vec![-1.0],
        vec![f64::NAN],
        vec![f64::INFINITY],
        vec![12.0; 4097],
    ] {
        assert!(gro_to_lammps_atomic(&text, &[1, 1], &masses, GroLammpsUnits::Real).is_err());
    }
    assert!(
        gro_to_lammps_atomic(
            &frame("2 3 4", [true, false]),
            &[1, 1],
            &[12.0],
            GroLammpsUnits::Real
        )
        .is_err()
    );
}

#[test]
fn empty_frame_and_large_input_limits_are_explicit() {
    let report =
        gro_to_lammps_atomic("empty\n0\n2 3 4\n", &[], &[12.0], GroLammpsUnits::Real).unwrap();
    assert!(report.data.atoms.is_empty());
    assert!(
        gro_to_lammps_atomic("too many\n100001\n", &[], &[12.0], GroLammpsUnits::Real).is_err()
    );
    assert!(
        gro_to_lammps_atomic(
            &"x".repeat((16 << 20) + 1),
            &[],
            &[12.0],
            GroLammpsUnits::Real
        )
        .is_err()
    );
}
