#![allow(clippy::redundant_closure)]

//! Per-descriptor timing over a SMILES file (first N rows).
use chematic_chem as c;
use std::time::Instant;
fn main() {
    let args: Vec<String> = std::env::args().collect();
    let n: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(200);
    let smis: Vec<String> = std::fs::read_to_string(&args[1])
        .unwrap()
        .lines()
        .filter_map(|l| l.split_whitespace().next().map(str::to_string))
        .take(n)
        .collect();
    macro_rules! t {
        ($name:expr, $f:expr) => {{
            let mols: Vec<_> = smis
                .iter()
                .filter_map(|s| chematic_smiles::parse(s).ok())
                .collect();
            let st = Instant::now();
            for m in &mols {
                std::hint::black_box($f(m));
            }
            println!(
                "{:28} {:10.1} us/mol",
                $name,
                st.elapsed().as_secs_f64() * 1e6 / mols.len() as f64
            );
        }};
    }
    t!("gate", |m| chematic_perception::rdkit_model_may_disagree(m));
    t!("h_suppressed", |m| {
        chematic_smiles::rdkit_hydrogen_suppressed(m).is_some()
    });
    t!("parity_arom", |m| {
        chematic_perception::with_rdkit_parity_view(m, |v| v.is_ok())
    });
    t!("parity_identity", |m| {
        chematic_perception::rdkit_parity_view_is_identity(m)
    });
    t!("sssr_count", |m| chematic_perception::sssr_ring_count(m));
    t!("ring_flags", |m| {
        chematic_perception::ring_bond_flags_shared(m).len()
    });
    t!("ring_bundle", |m| c::ring_bundle(m).ring_count);
    t!("logp_and_mr", |m| c::logp_and_mr(m));
    t!("tpsa", |m| c::tpsa(m));
    t!("rdkit_tpsa", |m| c::rdkit_tpsa(m));
    t!("hbd", |m| c::hbd_count(m));
    t!("num_stereocenters", |m| c::num_stereocenters(m));
    t!("fsp3", |m| c::fsp3(m));
    t!("qed", |m| c::qed(m));
    t!("sa_score", |m| c::sa_score(m));
    t!("labute_asa", |m| c::labute_asa(m));
    t!("bertz_ct", |m| c::bertz_ct(m));
    t!("wiener", |m| c::wiener_index(m));
    t!("schultz", |m| c::schultz_mti(m));
    t!("gutman", |m| c::gutman_mti(m));
    t!("vabc", |m| c::vabc(m));
    t!("grav", |m| c::gravitational_index(m));
    t!("kappa1", |m| c::kappa1(m));
    t!("kappa2", |m| c::kappa2(m));
    t!("kappa3", |m| c::kappa3(m));
    t!("chi_all", |m| c::chi_all(m));
    t!("chi1v", |m| c::chi1v(m));
    t!("num_unspec", |m| c::num_unspecified_stereocenters(m));
    t!("estate_all", |m| c::estate_all(m));
    t!("ghose", |m| c::ghose_passes(m));
    t!("reos", |m| c::reos_passes(m));
    t!("pains", |m| c::pains_passes(m));
    t!("ro3", |m| c::ro3_passes(m));
    t!("lead_like", |m| c::lead_like_passes(m));
    t!("pfizer", |m| c::pfizer_3_75_passes(m));
    t!("pka_base", |m| c::pka_base(m));
    t!("pka_acid", |m| c::pka_acid(m));
    t!("mcf", |m| c::mcf_passes(m));
    t!("slogp_vsa", |m| c::slogp_vsa(m));
    t!("smr_vsa", |m| c::smr_vsa(m));
    t!("peoe_vsa", |m| c::peoe_vsa(m));
    t!("estate_vsa", |m| c::estate_vsa(m));
    t!("num_rings", |m| c::rdkit_num_rings(m));
    t!("rdkit_sssr_full", |m| {
        chematic_perception::rdkit_sssr_ring_order(m).map(|r| r.len())
    });
    t!("om_dative", |m| {
        chematic_perception::rdkit_organometallic_dative_bonds(m).len()
    });
    t!("formula", |m| c::calc_mol_formula(m));
}
