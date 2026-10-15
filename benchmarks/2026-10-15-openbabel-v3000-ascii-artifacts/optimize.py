from pathlib import Path
import json,hashlib
s=Path('/private/tmp/chematic-openbabel-ascii-source-20261015');r=Path('/private/tmp/chematic-openbabel-ascii-20261015');p=s/'crates/chematic-mol/src/mol3000.rs';old=p.read_text();a=old.index('pub fn read_mol_v3000_with_diagnostics(');b=old.index('/// Parse a MOL V3000 string into',a);reference=old[a:b].replace('pub fn read_mol_v3000_with_diagnostics(', 'fn reference_read(',1)
needle='        tokens.extend(line.payload.split_whitespace());';assert old.count(needle)==1
new='''        // ASCII splitting has the same token boundaries except for vertical
        // tab, which Rust excludes from ASCII whitespace. Retain the Unicode
        // path for that character and for every non-ASCII payload.
        if line.payload.is_ascii() && !line.payload.contains('\\u{b}') {
            tokens.extend(line.payload.split_ascii_whitespace());
        } else {
            tokens.extend(line.payload.split_whitespace());
        }'''
t=old.replace(needle,new,1)+'\n#[cfg(test)]\n#[path = "mol3000/ascii_whitespace_tests.rs"]\nmod ascii_whitespace_tests;\n';p.write_text(t)
fixture='benchmarks/fixtures/ethanol.v3000';separators=list(range(128))+[0x85,0xa0,0x1680,*range(0x2000,0x200b),0x2028,0x2029,0x202f,0x205f,0x3000,0x200b,0xfeff]
assert len(separators)==len(set(separators))==149
recipe={'fixture':fixture,'separator_codepoints':separators,'framing':['direct','continued'],'document_variants':['original','invalid_atom_index','invalid_coordinate'],'expected_cases':894};payload=json.dumps(recipe,indent=2)+'\n';q=s/'validation/benchmark_corpora/v3000-ascii-whitespace-v1.recipe.json';q.write_text(payload)
manifest={'id':'v3000-ascii-whitespace-v1','path':q.relative_to(s).as_posix(),'sha256':hashlib.sha256(q.read_bytes()).hexdigest(),'provenance':{'source':'Original token-separator recipe over the existing project ethanol V3000 fixture; no new molecular records or third-party data','generated':'2026-10-15','source_path':fixture,'source_sha256':hashlib.sha256((s/fixture).read_bytes()).hexdigest(),'transformation':'Each separator replaces spaces after the exact V30 prefix. Repeat with continuations and malformed atom index/coordinate fields. Original headers and chemical records remain.'},'license':'MIT OR Apache-2.0','redistribution':{'status':'allowed','rationale':'Project-authored syntax recipe and existing project fixture'},'minimization':'All ASCII characters, all remaining Unicode whitespace, and two Unicode non-whitespace boundaries; one existing ethanol record suffices for tokenization and typed-error comparisons. Broad chemistry remains in the unchanged 2064-record corpus.','oracle':{'name':'Frozen original read_mol_v3000_with_diagnostics with split_whitespace','version':'source 1.1.4 plus eleven accepted source candidates','mol3000_sha256':hashlib.sha256(old.encode()).hexdigest()},'expected':{'cases':894,'status':'Exact typed errors or full indexed molecular fields, coordinate bits, metadata, geometry/stereo diagnostics and rewrite bytes equal original','absolute_tolerance':0}};(q.with_name('v3000-ascii-whitespace-v1.manifest.json')).write_text(json.dumps(manifest,indent=2)+'\n')
test='//! Frozen reader oracle; provenance and recipe in v3000-ascii-whitespace-v1.manifest.json.\nuse super::*;\n\n'+reference+'''
fn compare(input: &str) {
    let expected = reference_read(input);
    let actual = read_mol_v3000_with_diagnostics(input);
    match (expected, actual) {
        (Err(left), Err(right)) => assert_eq!(format!("{left:?}"), format!("{right:?}")),
        (Ok(left), Ok(right)) => {
            assert_eq!(left.metadata, right.metadata);
            assert_eq!(left.mol.atoms().map(|(idx, atom)| format!("{idx:?}:{atom:?}:{:?}:{:?}:{:?}:{:?}", left.mol.stereo_neighbor_order(idx), left.mol.r_group_label(idx), left.mol.neighbors(idx).collect::<Vec<_>>(), left.mol.atom_tag(idx))).collect::<Vec<_>>(), right.mol.atoms().map(|(idx, atom)| format!("{idx:?}:{atom:?}:{:?}:{:?}:{:?}:{:?}", right.mol.stereo_neighbor_order(idx), right.mol.r_group_label(idx), right.mol.neighbors(idx).collect::<Vec<_>>(), right.mol.atom_tag(idx))).collect::<Vec<_>>());
            assert_eq!(left.mol.bonds().map(|pair| format!("{pair:?}")).collect::<Vec<_>>(), right.mol.bonds().map(|pair| format!("{pair:?}")).collect::<Vec<_>>());
            assert_eq!(left.coords.iter().map(|(x,y)| [x.to_bits(),y.to_bits()]).collect::<Vec<_>>(), right.coords.iter().map(|(x,y)| [x.to_bits(),y.to_bits()]).collect::<Vec<_>>());
            assert_eq!(format!("{:?}:{:?}:{:?}:{:?}:{:?}:{:?}:{:?}:{:?}",left.stereo_diagnostics,left.ez_diagnostics,left.conformer,left.coordinate_dimension,left.geometry_rank,left.stereo3d_diagnostics,left.square_planar_diagnostics,left.mol.stereo_groups()),format!("{:?}:{:?}:{:?}:{:?}:{:?}:{:?}:{:?}:{:?}",right.stereo_diagnostics,right.ez_diagnostics,right.conformer,right.coordinate_dimension,right.geometry_rank,right.stereo3d_diagnostics,right.square_planar_diagnostics,right.mol.stereo_groups()));
            assert_eq!(write_mol_v3000(&left.mol,&left.metadata,&left.coords),write_mol_v3000(&right.mol,&right.metadata,&right.coords));
        }
        _ => panic!("reader result changed"),
    }
}

#[test]
fn ascii_unicode_and_continued_reader_match_original() {
    let original = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../benchmarks/fixtures/ethanol.v3000"));
    let recipe: serde_json::Value = serde_json::from_str(include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../validation/benchmark_corpora/v3000-ascii-whitespace-v1.recipe.json"))).unwrap();
    for codepoint in 0..=0x10ffff {
        if let Some(character) = char::from_u32(codepoint) && character.is_whitespace() {
            assert!(recipe["separator_codepoints"].as_array().unwrap().iter().any(|v| v.as_u64() == Some(codepoint as u64)));
        }
    }
    let mut cases = 0;
    let mut documents = vec![original.to_string()];
    for bad in ["bad", "not-a-float"] {
        let mut in_atoms = false;
        let mut changed = false;
        let mut document = String::new();
        for line in original.lines() {
            if line == "M  V30 BEGIN ATOM" { in_atoms = true; }
            if in_atoms && !changed && line.strip_prefix(V30_PREFIX).is_some_and(|p| p.as_bytes()[0].is_ascii_digit()) {
                let mut fields: Vec<_> = line[V30_PREFIX.len()..].split_whitespace().collect();
                fields[if bad == "bad" { 0 } else { 2 }] = bad;
                document.push_str(V30_PREFIX); document.push_str(&fields.join(" ")); changed = true;
            } else { document.push_str(line); }
            document.push('\\n');
        }
        assert!(changed); documents.push(document);
    }
    for document in documents {
        for codepoint in recipe["separator_codepoints"].as_array().unwrap() {
            let separator = char::from_u32(codepoint.as_u64().unwrap() as u32).unwrap().to_string();
            for continued in [false, true] {
                let mut input = String::new();
                for line in document.lines() {
                    if let Some(payload) = line.strip_prefix(V30_PREFIX) {
                        input.push_str(V30_PREFIX);
                        if continued && let Some((first,rest)) = payload.split_once(' ') {
                            input.push_str(first); input.push_str("-\\n"); input.push_str(V30_PREFIX); input.push_str(&rest.replace(' ', &separator));
                        } else { input.push_str(&payload.replace(' ', &separator)); }
                    } else { input.push_str(line); }
                    input.push('\\n');
                }
                compare(&input); cases += 1;
            }
        }
    }
    assert_eq!(cases,recipe["expected_cases"].as_u64().unwrap());
}
''';test=test.replace('left.mol.atom_tag(idx)','left.mol.atom_tag(idx)').replace('right.mol.atom_tag(idx)','right.mol.atom_tag(idx)');(s/'crates/chematic-mol/src/mol3000/ascii_whitespace_tests.rs').write_text(test)
p=s/'CHANGELOG.md';p.write_text(p.read_text().replace('### Performance\n','### Performance\n\n- Split ordinary ASCII V3000 payloads with the byte whitespace iterator,\n  retaining Unicode tokenization for non-ASCII text and vertical tabs.\n',1));print('ASCII lexer candidate and 894-case frozen-reader oracle written')
