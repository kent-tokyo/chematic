use chematic_3d::xyz::{XyzError, XyzRealError, parse_xyz, parse_xyz_real};
use serde_json::Value;

fn packet() -> Value {
    serde_json::from_str(include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../validation/xyz_numeric_read/cases.json"
    )))
    .unwrap()
}

#[test]
fn shared_real_tokens_and_every_coordinate_axis_obey_the_contract() {
    let packet = packet();
    for group in ["valid", "invalid"] {
        for case in packet[group].as_array().unwrap() {
            let token = case["token"].as_str().unwrap();
            let parsed = parse_xyz_real(token);
            if group == "valid" {
                let bits = u64::from_str_radix(case["bits"].as_str().unwrap(), 16).unwrap();
                assert_eq!(parsed.unwrap().to_bits(), bits, "{case}");
            } else {
                let reason = match case["reason"].as_str().unwrap() {
                    "underflow" => XyzRealError::Underflow,
                    "nonfinite" => XyzRealError::NonFinite,
                    _ => XyzRealError::Invalid,
                };
                assert_eq!(parsed.unwrap_err(), reason, "{case}");
            }
            for axis in 0..3 {
                let mut row = ["0", "0", "0"];
                row[axis] = token;
                let text = format!("1\nnumber\nC {}\n", row.join(" "));
                match parse_xyz(&text) {
                    Ok((mol, coords)) => {
                        assert_eq!(group, "valid", "{case}");
                        assert_eq!(mol.atom_count(), 1);
                        let p = coords.points[0];
                        assert_eq!([p.x, p.y, p.z][axis].to_bits(), parsed.unwrap().to_bits());
                    }
                    Err(error) => {
                        assert_eq!(group, "invalid", "{case}");
                        assert_eq!(error, XyzError::InvalidLine(3));
                    }
                }
            }
        }
    }
}

#[test]
fn finite_binary64_bit_patterns_remain_readable() {
    let mut state = 0x6a09_e667_f3bc_c909u64;
    for _ in 0..4096 {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        let value = f64::from_bits(state);
        if value.is_finite() {
            let text = format!("{value:.17e}");
            assert_eq!(parse_xyz_real(&text).unwrap().to_bits(), value.to_bits());
        }
    }
}
