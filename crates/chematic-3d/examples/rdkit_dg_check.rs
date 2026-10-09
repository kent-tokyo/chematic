//! Plain-DG check against RDKit: each input line is a JSON object
//! `{"n":..,"bounds":[n*n raw],"coords":[n*3] | null}` from RDKit's
//! `GetMoleculeBoundsMatrix` + `EmbedMolecule(useExpTorsionAnglePrefs=False,
//! useBasicKnowledge=False, boundsMat=..., randomSeed=42)` on molecules
//! without chiral or tetrahedral sets.
use chematic_3d::rdkit_embed::{BoundsMatrix, ChiralSet, EmbedArgs, EmbedParams, embed_points};
use std::io::BufRead;

fn main() {
    let path = std::env::args().nth(1).expect("input");
    let (mut ok, mut bad) = (0, 0);
    for line in std::io::BufReader::new(std::fs::File::open(path).unwrap()).lines() {
        let line = line.unwrap();
        let v: serde_json::Value = serde_json::from_str(&line).unwrap();
        let n = v["n"].as_u64().unwrap() as usize;
        let raw: Vec<f64> = v["bounds"]
            .as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_str().unwrap().parse::<f64>().unwrap())
            .collect();
        let mmat = BoundsMatrix::from_raw(n, raw);
        let sets = |k: &str| -> Vec<ChiralSet> {
            v[k].as_array()
                .unwrap()
                .iter()
                .map(|c| {
                    let idx: Vec<usize> = c[0]
                        .as_array()
                        .unwrap()
                        .iter()
                        .map(|x| x.as_u64().unwrap() as usize)
                        .collect();
                    ChiralSet {
                        idx: [idx[0], idx[1], idx[2], idx[3], idx[4]],
                        vol_lower: c[1].as_f64().unwrap(),
                        vol_upper: c[2].as_f64().unwrap(),
                        in_fused_small_rings: c[3].as_bool().unwrap(),
                    }
                })
                .collect()
        };
        let chiral = sets("chiral");
        let tet = sets("tet");
        let u = |x: &serde_json::Value| x.as_u64().unwrap() as usize;
        let dbe: Vec<(usize, usize, usize)> = v["dbe"]
            .as_array()
            .unwrap()
            .iter()
            .map(|t| (u(&t[0]), u(&t[1]), u(&t[2])))
            .collect();
        let sdb: Vec<([usize; 4], i32)> = v["sdb"]
            .as_array()
            .unwrap()
            .iter()
            .map(|t| {
                (
                    [u(&t[0][0]), u(&t[0][1]), u(&t[0][2]), u(&t[0][3])],
                    t[1].as_i64().unwrap() as i32,
                )
            })
            .collect();
        let args = EmbedArgs {
            mmat: &mmat,
            chiral_centers: &chiral,
            tetrahedral_centers: &tet,
            double_bond_ends: &dbe,
            stereo_double_bonds: &sdb,
            exp_torsions: None,
        };
        let got = embed_points(n, &args, &EmbedParams::default(), 42);
        let expect: Option<Vec<f64>> = v["coords"].as_array().map(|a| {
            a.iter()
                .map(|x| x.as_str().unwrap().parse::<f64>().unwrap())
                .collect()
        });
        let same = match (&got, &expect) {
            (Ok(Some(g)), Some(e)) => {
                g.iter()
                    .zip(e)
                    .map(|(a, b)| (a - b).abs())
                    .fold(0.0, f64::max)
                    <= 1e-8
            }
            (Ok(None), None) => true,
            _ => false,
        };
        if same {
            ok += 1;
        } else {
            bad += 1;
            if bad <= 100 {
                let md = match (&got, &expect) {
                    (Ok(Some(g)), Some(e)) => g
                        .iter()
                        .zip(e)
                        .map(|(a, b)| (a - b).abs())
                        .fold(0.0, f64::max),
                    _ => f64::NAN,
                };
                eprintln!(
                    "MISMATCH {} maxdiff={md:e} chiral={} tet={} got_some={:?}",
                    v["smiles"],
                    chiral.len(),
                    tet.len(),
                    got.as_ref().map(|x| x.is_some())
                );
            }
        }
    }
    println!("match {ok} mismatch {bad}");
}
