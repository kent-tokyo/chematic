from pathlib import Path
import json,hashlib
s=Path('/private/tmp/chematic-openbabel-zero-z-source-20261015');r=Path('/private/tmp/chematic-openbabel-zero-z-20261015');p=s/'crates/chematic-mol/src/mol3000.rs';old=p.read_text();a=old.index('fn parse_v3000_coords(');b=old.index('fn v3000_atom_identity(',a);reference=old[a:b].replace('fn parse_v3000_coords(', 'fn reference_coords(',1)
needle='    Ok((parse(2, "x")?, parse(3, "y")?, parse(4, "z")?))'
assert needle in old
t=old.replace(needle,'''    Ok((
        parse(2, "x")?,
        parse(3, "y")?,
        // Common finite positive zeros need no general float parsing.
        // Keep every other spelling, including negative zero, on that path.
        match tokens[4] {
            "0" | "0.0000" => 0.0,
            _ => parse(4, "z")?,
        },
    ))''',1)
t+='\n#[cfg(test)]\n#[path = "mol3000/zero_z_tests.rs"]\nmod zero_z_tests;\n';p.write_text(t)
recipe={'id':'v3000-zero-z-lexemes-v1','lexemes':['0','0.0000','-0','-0.0000','+0.0000','0.0','0e0','0e999','0E-999','.0','-0e-999','1','-1','1e-300','1e300','0.0001','NaN','inf','-inf','1e999','','bad',' 0','0 '],'lines':[0,1,42,18446744073709551615],'cartesian':'All x/y/z lexeme combinations at all listed source lines','random_f64':{'seed_hex':'5eed1234abcd9876','algorithm':'x ^= x << 13; x ^= x >> 7; x ^= x << 17; fixed u64 wrapping; reinterpret as f64, to_string for z','count':10000}}
relative='validation/benchmark_corpora/v3000-zero-z-speed-v1.recipe.json';payload=(json.dumps(recipe,indent=2,sort_keys=True)+'\n').encode();(s/relative).write_bytes(payload)
manifest={'id':'v3000-zero-z-speed-v1','path':relative,'provenance':{'source':'Original synthetic coordinate-lexeme recipe; no third-party or new molecular data.','retrieved_or_generated':'2026-10-15','transformation':recipe['cartesian']+'; also use the deterministic random-f64 recipe for the general-parser fallback.'},'license':'MIT OR Apache-2.0','redistribution':{'status':'allowed','rationale':'Original project-authored syntax cases and deterministic numeric recipe.'},'minimization':'24 lexemes cover positive/negative zero, equivalent alternate spellings, finite extremes, nonfinite values and malformed values. Their Cartesian product checks x/y/z failure precedence. Four source lines check retained error locations.','oracle':{'name':'Frozen pre-change chematic parse_v3000_coords','version':'1.1.4 plus seven accepted source candidates','settings':'Exact tuple float bits or complete typed-error Debug representation; no oracle refresh.','mol3000_sha256':hashlib.sha256((r/'baseline-mol3000.rs').read_bytes()).hexdigest()},'expected':{'coordinate_units':'angstrom','absolute_tolerance':0,'status':'All outputs and diagnostics exactly match the original implementation; negative-zero bits preserved.','cartesian_cases':24**3*4,'fallback_cases':10000},'sha256':hashlib.sha256(payload).hexdigest()}
(s/'validation/benchmark_corpora/v3000-zero-z-speed-v1.manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
test='''//! Differential oracle; see the v3000-zero-z fixture recipe and manifest.
use super::*;

'''+reference+'''
fn compare(x: &str, y: &str, z: &str, line: usize) {
    let tokens = ["1", "C", x, y, z];
    let expected = reference_coords(&tokens, line);
    let actual = parse_v3000_coords(&tokens, line);
    match (expected, actual) {
        (Ok((x,y,z)), Ok((a,b,c))) => assert_eq!([x.to_bits(),y.to_bits(),z.to_bits()], [a.to_bits(),b.to_bits(),c.to_bits()]),
        (Err(left), Err(right)) => assert_eq!(format!("{left:?}"), format!("{right:?}")),
        (left,right) => panic!("coordinate result differs: {left:?} vs {right:?}"),
    }
}

fn recipe() -> serde_json::Value {
    serde_json::from_str(include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../validation/benchmark_corpora/v3000-zero-z-speed-v1.recipe.json"))).unwrap()
}

#[test]
fn zero_spellings_and_failure_precedence_match_original() {
    let recipe = recipe();
    let values: Vec<_> = recipe["lexemes"].as_array().unwrap().iter().map(|v| v.as_str().unwrap()).collect();
    for line in recipe["lines"].as_array().unwrap() {
        // Bound the portable test location to the platform's usize range.
        let line = line.as_u64().unwrap().min(usize::MAX as u64) as usize;
        for &x in &values {
            for &y in &values {
                for &z in &values { compare(x,y,z,line); }
            }
        }
    }
    assert_eq!(values.len(),24);
}

#[test]
fn deterministic_f64_fallback_and_signed_zero_match_original() {
    let recipe = recipe();
    let random = &recipe["random_f64"];
    let mut bits = u64::from_str_radix(random["seed_hex"].as_str().unwrap(),16).unwrap();
    for _ in 0..random["count"].as_u64().unwrap() {
        bits ^= bits << 13;
        bits ^= bits >> 7;
        bits ^= bits << 17;
        let value = f64::from_bits(bits).to_string();
        compare("1", "-1", &value, 42);
    }
}
'''
(s/'crates/chematic-mol/src/mol3000/zero_z_tests.rs').write_text(test)
p=s/'CHANGELOG.md';t=p.read_text().replace('### Performance\n','### Performance\n\n- Parse the common V3000 z-coordinate spellings `0` and `0.0000` directly\n  as positive zero, retaining the general parser for every other spelling,\n  signed-zero bits and existing coordinate-error precedence.\n',1);p.write_text(t)
print('Positive-zero z fast path prepared')
