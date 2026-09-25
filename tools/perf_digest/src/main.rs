//! Output-digest and cold-timing harness for chematic performance work.
//!
//! * `digest <out.tsv> <corpus.smi>...` — one line per (corpus row, op) with the
//!   op's full output; every op runs on a fresh clone so no memoized state
//!   leaks between ops. Two revisions' digests must be byte-identical for an
//!   output-preserving change (see `scripts/perf_digest_diff.sh`).
//! * `time <corpus.smi>` — cold µs/mol per op (fresh clones every repeat;
//!   `REPS`, default 3, best-of), plus a shared-molecule session total.
//! * `patterns <smarts.txt> <corpus.smi>` — per-SMARTS existence-search cost.
//!
//! `ONLY=a,b,c` restricts the op set.
use chematic_core::Molecule;
use std::fmt::Write as _;
use std::time::Instant;

fn load(path: &str) -> Vec<Molecule> {
    std::fs::read_to_string(path)
        .unwrap()
        .lines()
        .filter_map(|l| l.split_whitespace().next())
        .filter_map(|s| chematic_smiles::parse(s).ok())
        .collect()
}

fn smarts_queries() -> &'static Vec<chematic_smarts::QueryMolecule> {
    static Q: std::sync::OnceLock<Vec<chematic_smarts::QueryMolecule>> = std::sync::OnceLock::new();
    Q.get_or_init(|| {
        [
            "c1ccccc1",
            "[OH]",
            "C(=O)N",
            "[#7;R]",
            "c1ccc2ccccc2c1",
            "[CX3](=O)[OX2H1]",
            "[NX3;H2,H1;!$(NC=O)]",
            "*~*~*~*~*~*",
            "[#6]~[#7]",
            "[R2]",
            "[r5]",
            "[$(C=O)]O",
            "[!#6;!#1]~*~[!#6;!#1]",
            "c:c-[CH3]",
            "[C;X4;!R]-[N;R]",
        ]
        .iter()
        .map(|s| chematic_smarts::parse_smarts(s).unwrap())
        .collect()
    })
}

type Op = (&'static str, fn(&Molecule) -> String);

fn ops() -> Vec<Op> {
    #[allow(unused_mut)]
    let mut ops: Vec<Op> = vec![
        ("sssr", |m| {
            format!("{:?}", chematic_perception::find_sssr(m).rings())
        }),
        ("ring_count", |m| chematic_chem::ring_count(m).to_string()),
        ("sssr_n", |m| {
            // SSSR cost without formatting the rings.
            let r = chematic_perception::find_sssr_shared(m);
            let sum: usize = r.rings().iter().map(|x| x.len()).sum();
            format!("{} {}", r.ring_count(), sum)
        }),
        ("aromatic_ring_count", |m| {
            chematic_chem::aromatic_ring_count(m).to_string()
        }),
        ("logp", |m| {
            format!("{:.12}", chematic_chem::logp_crippen(m))
        }),
        ("mr", |m| {
            format!("{:.12}", chematic_chem::molar_refractivity(m))
        }),
        ("tpsa", |m| format!("{:.12}", chematic_chem::tpsa(m))),
        ("hba", |m| chematic_chem::hba_count(m).to_string()),
        ("hbd", |m| chematic_chem::hbd_count(m).to_string()),
        ("rotb", |m| {
            chematic_chem::rotatable_bond_count(m).to_string()
        }),
        ("ring_bundle", |m| {
            let b = chematic_chem::ring_bundle(m);
            format!(
                "{} {} {} {}",
                b.ring_count, b.aromatic_ring_count, b.hba_count, b.rotatable_bond_count
            )
        }),
        ("qed", |m| format!("{:.12}", chematic_chem::qed(m))),
        ("ecfp4", |m| format!("{:?}", chematic_fp::ecfp4(m))),
        (
            "rdkit_ecfp4",
            |m| match chematic_fp::rdkit_morgan_ecfp4_experimental(m) {
                Ok(r) => format!("{:?}", r.fingerprint),
                Err(e) => format!("ERR {e}"),
            },
        ),
        ("canonical", |m| chematic_smiles::canonical_smiles(m)),
        ("stereocenters", |m| {
            chematic_chem::num_stereocenters(m).to_string()
        }),
        ("inchi", |m| chematic_inchi::inchi(m)),
        ("smarts_matches", |m| {
            let mut out = String::new();
            for q in smarts_queries() {
                let _ = write!(out, "{:?};", chematic_smarts::find_matches(q, m));
            }
            out
        }),
        ("smarts_nouniq", |m| {
            let cfg = chematic_smarts::MatchConfig {
                uniquify: false,
                ..Default::default()
            };
            let mut out = String::new();
            for q in smarts_queries() {
                let _ = write!(
                    out,
                    "{:?};",
                    chematic_smarts::find_matches_with_config(q, m, &cfg)
                );
            }
            out
        }),
        ("has_sub", |m| {
            smarts_queries()
                .iter()
                .map(|q| {
                    if chematic_smarts::has_match_bounded(
                        q,
                        m,
                        &chematic_perception::find_sssr(m),
                        &Default::default(),
                    ) == chematic_smarts::MatchOutcome::Found
                    {
                        '1'
                    } else {
                        '0'
                    }
                })
                .collect()
        }),
        ("has_sub_first", |m| {
            // First-embedding existence form used by the Python/bulk screens.
            let cfg = chematic_smarts::MatchConfig {
                max_matches: Some(1),
                uniquify: false,
                ..Default::default()
            };
            smarts_queries()
                .iter()
                .map(|q| {
                    if chematic_smarts::find_matches_with_config(q, m, &cfg).is_empty() {
                        '0'
                    } else {
                        '1'
                    }
                })
                .collect()
        }),
        ("rdkit_parity_view", |m| {
            match chematic_perception::apply_aromaticity_rdkit_parity_experimental(m) {
                Ok(p) => {
                    let atoms: String = p
                        .atoms()
                        .map(|(_, a)| if a.aromatic { 'a' } else { '.' })
                        .collect();
                    let bonds: String = p.bonds().map(|(_, b)| format!("{:?},", b.order)).collect();
                    format!("{atoms}|{bonds}")
                }
                Err(e) => format!("ERR {e}"),
            }
        }),
        ("rdkit_tpsa", |m| format!("{:.12}", chematic_chem::rdkit_tpsa(m))),
        ("rdkit_ecfp4_bitinfo", |m| match chematic_fp::rdkit_morgan_ecfp4_experimental(m) {
            Ok(r) => format!("{:?} {:?}", r.fingerprint, {
                let mut v: Vec<_> = r.raw_bit_info.iter().collect();
                v.sort();
                v
            }),
            Err(e) => format!("ERR {e}"),
        }),
        ("largest_frag", |m| {
            let f = chematic_chem::largest_fragment(m);
            let atoms: String = f
                .atoms()
                .map(|(_, a)| format!("{}{}{},", a.element.symbol(), a.charge, a.aromatic as u8))
                .collect();
            let bonds: String = f
                .bonds()
                .map(|(_, b)| format!("{}-{}:{:?},", b.atom1.0, b.atom2.0, b.order))
                .collect();
            format!(
                "{} | {} | {}",
                chematic_smiles::canonical_smiles(&f),
                atoms,
                bonds
            )
        }),
        ("standardize", |m| {
            let f = chematic_chem::standardize(m, &chematic_chem::StandardizeOptions::default());
            let atoms: String = f
                .atoms()
                .map(|(_, a)| {
                    format!(
                        "{}{}{}{}{:?},",
                        a.element.symbol(),
                        a.charge,
                        a.aromatic as u8,
                        a.isotope.unwrap_or(0),
                        a.chirality
                    )
                })
                .collect();
            let bonds: String = f
                .bonds()
                .map(|(_, b)| format!("{}-{}:{:?},", b.atom1.0, b.atom2.0, b.order))
                .collect();
            format!(
                "{} | {} | {}",
                chematic_smiles::canonical_smiles(&f),
                atoms,
                bonds
            )
        }),
        ("cip", |m| {
            format!("{:?}", chematic_chem::assign_cip(m).assignments)
        }),
        ("stereo_idx", |m| {
            format!("{:?}", chematic_chem::potential_stereocenter_indices(m))
        }),
        ("maccs", |m| format!("{:?}", chematic_fp::maccs::maccs(m))),
        ("pains", |m| chematic_chem::pains_passes(m).to_string()),
        ("brenk", |m| chematic_chem::brenk_passes(m).to_string()),
        // Binding-surface (Python/WASM) SMARTS semantics: perceived aromatic view.
        ("has_sub_perceived", |m| {
            let cfg = chematic_smarts::MatchConfig {
                max_matches: Some(1),
                uniquify: false,
                ..Default::default()
            };
            smarts_queries()
                .iter()
                .map(|q| {
                    if chematic_smarts::has_match_perceived(q, m, &cfg) {
                        '1'
                    } else {
                        '0'
                    }
                })
                .collect()
        }),
        ("has_sub_perceived_1", |m| {
            // One query per fresh molecule: the Python per-call shape.
            let cfg = chematic_smarts::MatchConfig {
                max_matches: Some(1),
                uniquify: false,
                ..Default::default()
            };
            chematic_smarts::has_match_perceived(&smarts_queries()[1], m, &cfg).to_string()
        }),
        ("find_matches_perceived", |m| {
            let mut out = String::new();
            for q in smarts_queries() {
                let mut v: Vec<Vec<u32>> = chematic_smarts::find_matches_perceived(
                    q,
                    m,
                    &chematic_smarts::MatchConfig::default(),
                )
                .into_iter()
                .map(|mm| {
                    let mut k: Vec<(usize, u32)> = mm.into_iter().map(|(a, b)| (a, b.0)).collect();
                    k.sort_unstable();
                    k.into_iter().map(|(_, b)| b).collect()
                })
                .collect();
                v.sort();
                let _ = write!(out, "{v:?};");
            }
            out
        }),
        ("match_atom_sets_perceived", |m| {
            // Python `Mol.find_matches` output. The candidate build computes it
            // with the candidate API; the base build from the maps.
            let mut out = String::new();
            for q in smarts_queries() {
                let cfg = chematic_smarts::MatchConfig::default();
                #[cfg(feature = "candidate-apis")]
                let v = chematic_smarts::find_match_atom_sets_perceived(q, m, &cfg);
                #[cfg(not(feature = "candidate-apis"))]
                let v = chematic_smarts::find_matches_perceived(q, m, &cfg)
                    .into_iter()
                    .map(|mm| {
                        let mut v: Vec<usize> = mm.values().map(|a| a.0 as usize).collect();
                        v.sort_unstable();
                        v
                    })
                    .collect::<Vec<_>>();
                let _ = write!(out, "{v:?};");
            }
            out
        }),
        ("rdkit_parity_full", |m| {
            // Every atom/bond field of the RDKit-parity view plus H counts and
            // the canonical SMILES written from it.
            match chematic_perception::apply_aromaticity_rdkit_parity_shared(m).as_ref() {
                Ok(p) => {
                    let mut out = String::new();
                    for (i, a) in p.atoms() {
                        let _ = write!(
                            out,
                            "{:?}/{}/{}/{:?}/{}/{:?}/{:?},",
                            a.element,
                            a.charge,
                            a.aromatic as u8,
                            a.hydrogen_count,
                            chematic_core::implicit_hcount(p, i),
                            a.chirality,
                            p.neighbors(i).collect::<Vec<_>>()
                        );
                    }
                    for (bi, b) in p.bonds() {
                        let _ = write!(
                            out,
                            "{}-{}:{:?}:{:?},",
                            b.atom1.0,
                            b.atom2.0,
                            b.order,
                            p.bond_direction(bi)
                        );
                    }
                    let _ = write!(out, "|{}", chematic_smiles::canonical_smiles(p));
                    out
                }
                Err(e) => format!("ERR {e}"),
            }
        }),
        ("rdkit_ecfp4_bits", |m| match chematic_fp::rdkit_morgan_ecfp4_bitvec(m) {
            Ok(f) => format!("{f:?}"),
            Err(e) => format!("ERR {e}"),
        }),
        ("rdkit_ecfp4_prepared", |m| match chematic_fp::prepare_rdkit_morgan_ecfp4(m) {
            Ok(p) => format!("{:?}", p.bitvec()),
            Err(e) => format!("ERR {e}"),
        }),
        ("kekulize", |m| match chematic_core::kekulize(m) {
            Ok(k) => {
                let mut v: Vec<_> = k.into_iter().map(|(b, o)| (b.0, o)).collect();
                v.sort_by_key(|&(b, _)| b);
                format!("{v:?}")
            }
            Err(e) => format!("ERR {e}"),
        }),
        ("rdkit_parity_ok", |m| {
            match chematic_perception::apply_aromaticity_rdkit_parity_shared(m).as_ref() {
                Ok(p) => p.atom_count().to_string(),
                Err(e) => format!("ERR {e}"),
            }
        }),
        ("rdkit_rdk_fp", |m| format!("{:?}", chematic_fp::rdkit_rdk_fp(m))),
        ("rdkit_pattern_fp", |m| format!("{:?}", chematic_fp::rdkit_pattern_fp(m))),
        ("rdkit_atom_pair", |m| format!("{:?}", chematic_fp::rdkit_atom_pair_fp(m))),
        ("rdkit_torsion", |m| format!("{:?}", chematic_fp::rdkit_torsion_fp(m))),
        ("chi_each", |m| {
            format!(
                "{:?}",
                [
                    chematic_chem::chi0(m),
                    chematic_chem::chi1(m),
                    chematic_chem::chi2(m),
                    chematic_chem::chi3(m),
                    chematic_chem::chi4(m),
                    chematic_chem::chi0v(m),
                    chematic_chem::chi1v(m),
                    chematic_chem::chi2v(m),
                    chematic_chem::chi3v(m),
                    chematic_chem::chi4v(m),
                ]
            )
        }),
        ("chi1v", |m| format!("{:?}", chematic_chem::chi1v(m))),
        ("chi_all", |m| format!("{:?}", chematic_chem::chi_all(m))),
        ("kappa", |m| {
            format!(
                "{:?}",
                (
                    chematic_chem::kappa1(m),
                    chematic_chem::kappa2(m),
                    chematic_chem::kappa3(m)
                )
            )
        }),
        ("rdkit_aromatic_ring_count", |m| {
            chematic_chem::rdkit_aromatic_ring_count(m).to_string()
        }),
    ];
    #[cfg(feature = "reference-oracles")]
    ops.push(("sssr_ref", |m| {
        format!(
            "{:?}",
            chematic_perception::sssr::find_sssr_horton_reference(m).rings()
        )
    }));
    ops
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mode = args[1].as_str();
    let only = std::env::var("ONLY").ok();
    let sel: Vec<Op> = ops()
        .into_iter()
        .filter(|(n, _)| only.as_ref().is_none_or(|o| o.split(',').any(|x| x == *n)))
        .collect();
    match mode {
        "time" => {
            let mols = load(&args[2]);
            let reps: usize = std::env::var("REPS")
                .ok()
                .and_then(|s| s.parse().ok())
                .unwrap_or(3);
            for (name, f) in &sel {
                let mut best = f64::MAX;
                for _ in 0..reps {
                    let fresh = mols.clone(); // cold: empty derived caches
                    let t = Instant::now();
                    for m in &fresh {
                        std::hint::black_box(f(m));
                    }
                    best = best.min(t.elapsed().as_secs_f64());
                }
                println!(
                    "{:22} {:10.2} µs/mol (cold)",
                    name,
                    best / mols.len() as f64 * 1e6
                );
            }
            if sel.len() > 1 {
                let fresh = mols.clone();
                let t = Instant::now();
                for m in &fresh {
                    for (_, f) in &sel {
                        std::hint::black_box(f(m));
                    }
                }
                println!(
                    "{:22} {:10.2} µs/mol (all selected ops, shared mol)",
                    "SESSION",
                    t.elapsed().as_secs_f64() / mols.len() as f64 * 1e6
                );
            }
        }
        "digest" => {
            // Each op computed on a fresh clone so no cache leaks between ops.
            // Streamed to disk: full SMARTS match maps make digests large.
            use std::io::Write as _;
            let file = std::fs::File::create(&args[2]).unwrap();
            let mut out = std::io::BufWriter::new(file);
            for path in &args[3..] {
                let mols = load(path);
                for (i, m) in mols.iter().enumerate() {
                    for (name, f) in &sel {
                        let fresh = m.clone();
                        writeln!(out, "{path}\t{i}\t{name}\t{}", f(&fresh)).unwrap();
                    }
                }
            }
            out.flush().unwrap();
        }
        "patterns" => {
            let mols = load(&args[3]);
            let cfg = chematic_smarts::MatchConfig {
                max_matches: Some(1),
                ..Default::default()
            };
            let mut rows = vec![];
            for p in std::fs::read_to_string(&args[2]).unwrap().lines() {
                let Ok(q) = chematic_smarts::parse_smarts(p) else {
                    continue;
                };
                let rings: Vec<_> = mols
                    .iter()
                    .map(|m| chematic_perception::find_sssr(m))
                    .collect();
                let t = Instant::now();
                let mut hits = 0;
                for (m, r) in mols.iter().zip(&rings) {
                    if !chematic_smarts::find_matches_with_rings_and_config(&q, m, r, &cfg)
                        .is_empty()
                    {
                        hits += 1;
                    }
                }
                rows.push((
                    t.elapsed().as_secs_f64() / mols.len() as f64 * 1e6,
                    p.to_string(),
                    hits,
                ));
            }
            rows.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
            let total: f64 = rows.iter().map(|r| r.0).sum();
            println!("total {total:.1} µs/mol over {} patterns", rows.len());
            for (t, p, h) in rows.iter().take(25) {
                println!("{t:8.2} µs  hits={h:5}  {p}");
            }
        }
        _ => panic!("mode"),
    }
}
