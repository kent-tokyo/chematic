from pathlib import Path
import hashlib,json,re
s=Path('/private/tmp/chematic-openbabel-flat-geometry-source-20261015');r=Path('/private/tmp/chematic-openbabel-flat-geometry-20261015')
p=s/'crates/chematic-mol/src/mol2000.rs';old=p.read_text();start=old.index('pub(crate) fn classify_geometry_rank(');end=old.index('/// Parse the dimensional-code',start);frozen=old[start:end].replace('pub(crate) fn classify_geometry_rank(', 'fn reference_geometry_rank(',1)
constants='\n'.join(re.findall(r'^const (?:Z_EPS|COPLANAR_EPS|COLLINEAR_EPS): f64 = [^;]+;',old,re.M))
needle='''    if points.is_empty() {
        return GeometryRank::Indeterminate;
    }
    if points.iter().all(|p| p.z.abs() < Z_EPS) {
        return GeometryRank::FlatZero;
    }
'''
assert needle in old
t=old.replace(needle,'''    if let Some(rank) = flat_geometry_rank(points.iter().map(|point| point.z)) {
        return rank;
    }
''',1)
helper='''/// Resolve the empty/flat prefix without materializing 3D points. The
/// strict near-zero comparison also preserves the non-finite fallback.
pub(crate) fn flat_geometry_rank(z: impl Iterator<Item = f64>) -> Option<GeometryRank> {
    let mut z = z.peekable();
    if z.peek().is_none() {
        Some(GeometryRank::Indeterminate)
    } else if z.all(|value| value.abs() < Z_EPS) {
        Some(GeometryRank::FlatZero)
    } else {
        None
    }
}

'''
pos=t.index('/// Classify a point cloud');t=t[:pos]+helper+t[pos:];p.write_text(t)
p=s/'crates/chematic-mol/src/mol3000.rs';old=p.read_text();a=old.index('fn v3000_geometry(');b=old.index('/// Apply supported CTfile',a);ref=old[a:b].replace('fn v3000_geometry(', 'fn reference_v3000_geometry(',1).replace('crate::mol2000::classify_geometry_rank(&points)','reference_geometry_rank(&points)')
needle='''    let points = coords
        .iter()
        .zip(raw_z)
        .map(|(&(x, y), &z)| Point3::new(x, y, z))
        .collect::<Vec<_>>();
    let rank = crate::mol2000::classify_geometry_rank(&points);
'''
assert needle in old
t=old.replace(needle,'''    // `zip` would truncate to the shorter coordinate channel. Test the
    // same z values first; empty and near-zero clouds never keep points.
    let flat_rank = crate::mol2000::flat_geometry_rank(raw_z.iter().copied().take(coords.len()));
    let points = if flat_rank.is_some() {
        Vec::new()
    } else {
        coords.iter().zip(raw_z)
            .map(|(&(x, y), &z)| Point3::new(x, y, z)).collect()
    };
    let rank = flat_rank.unwrap_or_else(|| crate::mol2000::classify_geometry_rank(&points));
''',1)
t+='\n#[cfg(test)]\n#[path = "mol3000/flat_geometry_tests.rs"]\nmod flat_geometry_tests;\n';p.write_text(t)
recipe={'id':'flat-geometry-prefix-v1','xy':'point i = (i, i mod 2)','z_pattern':'value[(i + offset) mod values.len()]','lengths':[0,1,2,3,4,5,6],'dimensions':['Unknown','TwoD','ThreeD'],'z_f64_bits_hex':['0000000000000000','8000000000000000','3f1a36e2eb1c432c','3f1a36e2eb1c432d','3f1a36e2eb1c432e','bf1a36e2eb1c432c','bf1a36e2eb1c432d','bf1a36e2eb1c432e','3ff0000000000000','bff0000000000000','7ff0000000000000','fff0000000000000','7ff8000000000001']}
# Boundary values are previous/exact/next representable f64 around 1e-4.
import struct,math
assert recipe['z_f64_bits_hex'][3]==struct.pack('>d',1e-4).hex()
relative='validation/benchmark_corpora/flat-geometry-speed-v1.recipe.json';payload=(json.dumps(recipe,indent=2,sort_keys=True)+'\n').encode();(s/relative).write_bytes(payload)
manifest={'id':'flat-geometry-speed-v1','path':relative,'provenance':{'source':'Original deterministic numeric recipe for differential testing, no new molecular records. Molecular comparisons reuse all rows of rdkit-js-browser-10k-v1.smi.','retrieved_or_generated':'2026-10-15','transformation':'Cross every listed coordinate-channel length, z offset and declared dimension; compare exact old classifier and V3000 geometry outputs. Also use full-length generated coordinates on every existing corpus molecule.'},'license':'MIT OR Apache-2.0','redistribution':{'status':'allowed','rationale':'Original project-authored recipe, no third-party numeric data; existing molecular corpus attribution and license remain unchanged.'},'minimization':'Lengths 0 through 6 cover empty, truncated, short, collinear and plane-defining inputs; z values cover signed zero, strict threshold neighbors and nonfinite fallback. Full original corpus is retained without filtering.','oracle':{'name':'Frozen pre-change chematic source candidate','version':'1.1.4 plus seven previously accepted candidates','settings':'Exact original classify_geometry_rank and v3000_geometry; source hashes identify the oracle, not current main.','mol2000_sha256':hashlib.sha256((r/'baseline-mol2000.rs').read_bytes()).hexdigest(),'mol3000_sha256':hashlib.sha256((r/'baseline-mol3000.rs').read_bytes()).hexdigest()},'expected':{'status':'Exact ranks, diagnostics, optional conformer point bits and graph annotations match frozen original; no tolerance change','coordinate_units':'angstrom','absolute_tolerance':0},'sha256':hashlib.sha256(payload).hexdigest()}
(s/'validation/benchmark_corpora/flat-geometry-speed-v1.manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
test='''//! Exact differential oracle and recipe provenance: see the flat-geometry
//! manifest under validation/benchmark_corpora.
use super::*;

'''+constants+'\n\n'+frozen+ref+'''
fn z_values() -> Vec<f64> {
    let recipe: serde_json::Value = serde_json::from_str(include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../validation/benchmark_corpora/flat-geometry-speed-v1.recipe.json"))).unwrap();
    recipe["z_f64_bits_hex"].as_array().unwrap().iter().map(|value| {
        f64::from_bits(u64::from_str_radix(value.as_str().unwrap(), 16).unwrap())
    }).collect()
}

fn compare(input: &Molecule, coords: &[(f64, f64)], z: &[f64], dimension: CoordinateDimension) {
    let mut expected_mol = input.clone();
    let mut actual_mol = input.clone();
    let expected = reference_v3000_geometry(&mut expected_mol, coords, z, dimension);
    let actual = v3000_geometry(&mut actual_mol, coords, z, dimension);
    assert_eq!(format!("{expected:?}"), format!("{actual:?}"));
    let bits = |conformer: &Option<Coords3D>| conformer.as_ref().map(|c| c.points.iter().map(|p| [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()]).collect::<Vec<_>>());
    assert_eq!(bits(&expected.0), bits(&actual.0));
    assert_eq!(expected_mol.atoms().collect::<Vec<_>>(), actual_mol.atoms().collect::<Vec<_>>());
    assert_eq!(expected_mol.bonds().collect::<Vec<_>>(), actual_mol.bonds().collect::<Vec<_>>());
    let points: Vec<_> = coords.iter().zip(z).map(|(&(x,y), &z)| Point3::new(x,y,z)).collect();
    assert_eq!(reference_geometry_rank(&points), crate::mol2000::classify_geometry_rank(&points));
}

#[test]
fn empty_truncated_threshold_and_nonfinite_channels_match_original() {
    let empty = MoleculeBuilder::new().build();
    let values = z_values();
    for xy_count in 0..=6 {
        let coords: Vec<_> = (0..xy_count).map(|i| (i as f64, (i % 2) as f64)).collect();
        for z_count in 0..=6 {
            for offset in 0..values.len() {
                let z: Vec<_> = (0..z_count).map(|i| values[(i+offset)%values.len()]).collect();
                for dimension in [CoordinateDimension::Unknown, CoordinateDimension::TwoD, CoordinateDimension::ThreeD] {
                    compare(&empty, &coords, &z, dimension);
                }
            }
        }
    }
}

#[test]
fn fixed_ten_thousand_graph_annotations_and_geometry_match_original() {
    let corpus = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi"));
    let values = z_values();
    let mut rows = 0;
    for line in corpus.lines().filter(|line| !line.trim().is_empty()) {
        let mol = chematic_smiles::parse_template(line.split_whitespace().next().unwrap()).unwrap();
        let coords: Vec<_> = (0..mol.atom_count()).map(|i| (i as f64, (i%2) as f64)).collect();
        for value in &values[..10] {
            let z = vec![*value; mol.atom_count()];
            for dimension in [CoordinateDimension::Unknown, CoordinateDimension::TwoD, CoordinateDimension::ThreeD] {
                compare(&mol, &coords, &z, dimension);
            }
        }
        rows += 1;
    }
    assert_eq!(rows, 10000);
}
'''
(s/'crates/chematic-mol/src/mol3000/flat_geometry_tests.rs').write_text(test)
p=s/'CHANGELOG.md';t=p.read_text().replace('### Performance\n','### Performance\n\n- Avoid allocating a temporary 3D point array when V3000 coordinates are\n  empty or flat near zero, reusing the existing strict geometry threshold\n  and preserving dimensional diagnostics and non-flat perception.\n',1);p.write_text(t)
print('Flat-geometry candidate prepared')
