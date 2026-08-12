//! Dumps `descriptors.rs` values that `scripts/descriptor_census.py` cannot
//! cheaply reach from Python for two different reasons:
//!
//! 1. No Python/WASM/MCP binding at all (confirmed by grepping chematic-py/
//!    chematic-mcp/chematic-wasm for their names — see
//!    docs/rfcs/descriptor_census_rfc.md): `moran_autocorr`, `geary_autocorr`,
//!    `information_content`, `mde_carbon`, and the plain `mmff94_charges`
//!    (shadowed in the Python API by `mmff94_charges_bci`).
//! 2. Only reachable via the monolithic `Mol.descriptors()` dict, which
//!    unconditionally computes ~130 other out-of-scope values too (QED, SA
//!    score, PAINS, drug_score, ...) — and on one pathological symmetric
//!    macrocycle in this census's corpus, `drug_score`'s PAINS/VF2
//!    substructure match takes several minutes (see docs/rfcs/descriptor_census_rfc.md's
//!    "VF2 performance" finding). `bcut2d` and `carbon_types` have no
//!    individual Python getter, so they are dumped here instead of paying
//!    that cost on every molecule.
//!
//! Writes one JSON object per line (JSONL) to stdout.
//!
//! Usage:
//!   cargo run -p chematic-chem --release --example descriptor_census_unbound \
//!       < scripts/descriptor_census_corpus.smi > validation/results/descriptor_census_unbound.jsonl

use chematic_chem::descriptors::{
    bcut2d, carbon_types, geary_autocorr, information_content, mde_carbon, mmff94_charges,
    moran_autocorr,
};
use chematic_core::Molecule;
use std::io::{self, BufRead, Write};
use std::panic;

fn fmt_vec(v: &[f64]) -> String {
    let parts: Vec<String> = v.iter().map(|x| format!("{x:.10}")).collect();
    format!("[{}]", parts.join(","))
}

/// Calls all 5 unbound functions, catching panics per-function so one
/// crashing descriptor doesn't hide the others' results for the same
/// molecule (diagnosis-only: a real panic found in `moran_autocorr` on this
/// corpus, see docs/rfcs/descriptor_census_rfc.md).
fn safe_call<T>(f: impl FnOnce() -> T + panic::UnwindSafe) -> Result<T, String> {
    panic::catch_unwind(f).map_err(|e| {
        if let Some(s) = e.downcast_ref::<&str>() {
            s.to_string()
        } else if let Some(s) = e.downcast_ref::<String>() {
            s.clone()
        } else {
            "panic (unknown payload)".to_string()
        }
    })
}

fn main() {
    // Suppress the default panic backtrace noise on stderr; we report panics
    // ourselves in the JSONL output instead.
    panic::set_hook(Box::new(|_| {}));

    let stdin = io::stdin();
    let stdout = io::stdout();
    let mut out = stdout.lock();
    for line in stdin.lock().lines() {
        let smi = line.unwrap();
        let smi = smi.trim();
        if smi.is_empty() {
            continue;
        }
        match chematic_smiles::parse(smi) {
            Ok(mol) => {
                let mol: Molecule = mol;
                let moran = safe_call(panic::AssertUnwindSafe(|| moran_autocorr(&mol)));
                let geary = safe_call(panic::AssertUnwindSafe(|| geary_autocorr(&mol)));
                let ic = safe_call(panic::AssertUnwindSafe(|| information_content(&mol)));
                let mde = safe_call(panic::AssertUnwindSafe(|| mde_carbon(&mol)));
                let mmff = safe_call(panic::AssertUnwindSafe(|| mmff94_charges(&mol)));
                let bc = bcut2d(&mol);
                let ct = carbon_types(&mol);

                let moran_json = match &moran {
                    Ok(v) => fmt_vec(v),
                    Err(e) => format!("null,\"moran_autocorr_error\":{:?}", e),
                };
                let geary_json = match &geary {
                    Ok(v) => fmt_vec(v),
                    Err(e) => format!("null,\"geary_autocorr_error\":{:?}", e),
                };
                let ic_json = match &ic {
                    Ok(v) => format!(
                        "{:.10},\"tic\":{:.10},\"sic\":{:.10},\"bic\":{:.10},\"cic\":{:.10}",
                        v.ic, v.tic, v.sic, v.bic, v.cic
                    ),
                    Err(e) => format!(
                        "null,\"tic\":null,\"sic\":null,\"bic\":null,\"cic\":null,\"ic_error\":{:?}",
                        e
                    ),
                };
                let mde_json = match &mde {
                    Ok(v) => fmt_vec(v),
                    Err(e) => format!("null,\"mde_carbon_error\":{:?}", e),
                };
                let mmff_json = match &mmff {
                    Ok(v) => {
                        let sum: f64 = v.iter().sum();
                        format!("{:.10},\"mmff94_charges_n\":{}", sum, v.len())
                    }
                    Err(e) => format!(
                        "null,\"mmff94_charges_n\":null,\"mmff94_charges_error\":{:?}",
                        e
                    ),
                };

                writeln!(
                    out,
                    "{{\"smiles\":{:?},\"parse_ok\":true,\"moran_autocorr\":{},\"geary_autocorr\":{},\"ic\":{},\"mde_carbon\":{},\"mmff94_charges_sum\":{},\
                     \"bcut2d_mwhi\":{:.10},\"bcut2d_mwlo\":{:.10},\"bcut2d_chghi\":{:.10},\"bcut2d_chglo\":{:.10},\
                     \"bcut2d_logphi\":{:.10},\"bcut2d_logplo\":{:.10},\"bcut2d_mrhi\":{:.10},\"bcut2d_mrlo\":{:.10},\
                     \"c1sp1\":{},\"c2sp1\":{},\"c1sp2\":{},\"c2sp2\":{},\"c3sp2\":{},\"c1sp3\":{},\"c2sp3\":{},\"c3sp3\":{}}}",
                    smi, moran_json, geary_json, ic_json, mde_json, mmff_json,
                    bc.mwhi, bc.mwlo, bc.chghi, bc.chglo, bc.logphi, bc.logplo, bc.mrhi, bc.mrlo,
                    ct.c1sp1, ct.c2sp1, ct.c1sp2, ct.c2sp2, ct.c3sp2, ct.c1sp3, ct.c2sp3, ct.c3sp3,
                )
                .unwrap();
            }
            Err(e) => {
                writeln!(
                    out,
                    "{{\"smiles\":{:?},\"parse_ok\":false,\"error\":{:?}}}",
                    smi,
                    e.to_string()
                )
                .unwrap();
            }
        }
    }
}
