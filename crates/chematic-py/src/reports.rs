//! `Report`/`compare` molecule-report bindings.

use crate::Mol;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;

/// A rendered HTML report returned by `chematic.report()` and `chematic.compare()`.
///
/// In Jupyter, writing ``report`` in a cell renders the HTML automatically.
/// Use ``report.save("path.html")`` to write to disk, or ``str(report)`` to get the HTML string.
///
///     report = chematic.report(mols, names=["aspirin", "ibuprofen"])
///     report.save("report.html")   # write to disk
///     display(report)              # Jupyter: renders inline
#[pyclass]
struct Report {
    html: String,
}

#[pymethods]
impl Report {
    fn _repr_html_(&self) -> &str {
        &self.html
    }

    fn save(&self, path: &str) -> PyResult<()> {
        std::fs::write(path, &self.html)
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))
    }

    fn __str__(&self) -> &str {
        &self.html
    }

    fn __repr__(&self) -> String {
        format!("Report({} bytes)", self.html.len())
    }
}

fn mol_report_to_dict<'py>(
    py: Python<'py>,
    report: &chematic_chem::MoleculeReport,
) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    d.set_item("canonical_smiles", &report.canonical_smiles)?;
    d.set_item("formula", &report.formula)?;
    match &report.murcko_scaffold_smiles {
        Some(s) => d.set_item("murcko_scaffold", s)?,
        None => d.set_item("murcko_scaffold", py.None())?,
    }
    let desc = PyDict::new(py);
    desc.set_item("mw", report.descriptors.molecular_weight)?;
    desc.set_item("exact_mass", report.descriptors.exact_mass)?;
    desc.set_item("tpsa", report.descriptors.tpsa)?;
    desc.set_item("logp", report.descriptors.logp)?;
    desc.set_item("molar_refractivity", report.descriptors.molar_refractivity)?;
    desc.set_item("hbd", report.descriptors.hbd)?;
    desc.set_item("hba", report.descriptors.hba)?;
    desc.set_item("rotatable_bonds", report.descriptors.rotatable_bonds)?;
    desc.set_item("heavy_atoms", report.descriptors.heavy_atom_count)?;
    desc.set_item("ring_count", report.descriptors.ring_count)?;
    desc.set_item("num_heteroatoms", report.descriptors.num_heteroatoms)?;
    desc.set_item("num_stereocenters", report.descriptors.num_stereocenters)?;
    desc.set_item("fsp3", report.descriptors.fsp3)?;
    desc.set_item("qed", report.descriptors.qed)?;
    desc.set_item("sa_score", report.descriptors.sa_score)?;
    desc.set_item("formal_charge", report.descriptors.formal_charge_sum)?;
    desc.set_item("labute_asa", report.descriptors.labute_asa)?;
    desc.set_item("bertz_ct", report.descriptors.bertz_ct)?;
    desc.set_item("wiener_index", report.descriptors.wiener_index)?;
    d.set_item("descriptors", desc)?;
    let filters = PyDict::new(py);
    filters.set_item("lipinski_passes", report.filters.lipinski_passes)?;
    filters.set_item("veber_passes", report.filters.veber_passes)?;
    filters.set_item("egan_passes", report.filters.egan_passes)?;
    filters.set_item("ghose_passes", report.filters.ghose_passes)?;
    filters.set_item("reos_passes", report.filters.reos_passes)?;
    filters.set_item("pains_passes", report.filters.pains_passes)?;
    let alerts: Vec<&str> = report
        .filters
        .pains_alerts
        .iter()
        .map(|s| s.as_str())
        .collect();
    filters.set_item("pains_alerts", alerts)?;
    d.set_item("filters", filters)?;
    let fgs: Vec<Bound<'py, PyDict>> = report
        .functional_groups
        .iter()
        .map(|fg| {
            let fd = PyDict::new(py);
            fd.set_item("name", &fg.name)?;
            fd.set_item("atom_indices", &fg.atom_indices)?;
            Ok(fd)
        })
        .collect::<PyResult<Vec<_>>>()?;
    d.set_item("functional_groups", fgs)?;
    let ngs: Vec<Bound<'py, PyDict>> = report
        .named_groups
        .iter()
        .map(|ng| {
            let nd = PyDict::new(py);
            nd.set_item("name", &ng.name)?;
            nd.set_item("atom_indices", &ng.atom_indices)?;
            Ok(nd)
        })
        .collect::<PyResult<Vec<_>>>()?;
    d.set_item("named_groups", ngs)?;
    Ok(d)
}

// ---------------------------------------------------------------------------
// Module-level functions
// ---------------------------------------------------------------------------

/// Generate a complete single-molecule analysis report.
///
/// Returns a dict with keys:
///   ``canonical_smiles``, ``formula``, ``murcko_scaffold``,
///   ``descriptors`` (dict), ``filters`` (dict),
///   ``functional_groups`` (list of dicts), ``named_groups`` (list of dicts).
///
///     report = chematic.molecule_report("CC(=O)Oc1ccccc1C(=O)O")
///     print(report["descriptors"]["mw"])
///     print(report["filters"]["lipinski_passes"])
#[pyfunction]
fn molecule_report<'py>(smiles: &str, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
    let report =
        chematic_chem::molecule_report(smiles).map_err(|e| PyValueError::new_err(e.to_string()))?;
    mol_report_to_dict(py, &report)
}

/// Screen a list of SMILES and return a batch report with diversity analysis.
///
/// Returns a dict with keys:
///   ``records`` (list of per-molecule dicts with ``input_index``, ``smiles``, ``report``, ``error``),
///   ``maxmin_picks`` (list of indices — most diverse subset),
///   ``butina_clusters`` (list of cluster index lists).
///
///     result = chematic.screen_smiles(smiles_list)
///     for rec in result["records"]:
///         if rec["error"] is None:
///             print(rec["report"]["descriptors"]["mw"])
#[pyfunction]
fn screen_smiles<'py>(smiles: Vec<String>, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
    let refs: Vec<&str> = smiles.iter().map(|s| s.as_str()).collect();
    let report = chematic_chem::screen_smiles(&refs);
    let d = PyDict::new(py);
    let records: Vec<Bound<'py, PyDict>> = report
        .records
        .iter()
        .map(|rec| {
            let r = PyDict::new(py);
            r.set_item("input_index", rec.input_index)?;
            r.set_item("smiles", &rec.input_smiles)?;
            match &rec.report {
                Some(mr) => {
                    r.set_item("report", mol_report_to_dict(py, mr)?)?;
                    r.set_item("error", py.None())?;
                }
                None => {
                    r.set_item("report", py.None())?;
                    r.set_item("error", rec.error.as_deref().unwrap_or("unknown error"))?;
                }
            }
            Ok(r)
        })
        .collect::<PyResult<Vec<_>>>()?;
    d.set_item("records", records)?;
    d.set_item("maxmin_picks", report.maxmin_picks)?;
    d.set_item("butina_clusters", report.butina_clusters)?;
    Ok(d)
}

/// Compare two or more SMILES strings and return pairwise similarity + descriptor deltas.
///
/// Returns a dict with keys:
///   ``reports`` (list of molecule report dicts),
///   ``pairwise`` (list of pairwise similarity dicts),
///   ``descriptor_deltas`` (list of delta dicts),
///   ``mcs_smiles`` (str or None).
///
///     result = chematic.compare_molecules(["c1ccccc1", "Cc1ccccc1"])
///     print(result["pairwise"][0]["ecfp4_tanimoto"])
#[pyfunction]
fn compare_molecules<'py>(smiles: Vec<String>, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
    let refs: Vec<&str> = smiles.iter().map(|s| s.as_str()).collect();
    let cmp = chematic_chem::compare_molecules(&refs)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let d = PyDict::new(py);
    let reports: Vec<Bound<'py, PyDict>> = cmp
        .reports
        .iter()
        .map(|r| mol_report_to_dict(py, r))
        .collect::<PyResult<Vec<_>>>()?;
    d.set_item("reports", reports)?;
    let pairwise: Vec<Bound<'py, PyDict>> = cmp
        .pairwise
        .iter()
        .map(|p| {
            let pd = PyDict::new(py);
            pd.set_item("left_index", p.left_index)?;
            pd.set_item("right_index", p.right_index)?;
            pd.set_item("ecfp4_tanimoto", p.similarities.ecfp4_tanimoto)?;
            pd.set_item("maccs_tanimoto", p.similarities.maccs_tanimoto)?;
            pd.set_item("atom_pair_tanimoto", p.similarities.atom_pair_tanimoto)?;
            Ok(pd)
        })
        .collect::<PyResult<Vec<_>>>()?;
    d.set_item("pairwise", pairwise)?;
    let deltas: Vec<Bound<'py, PyDict>> = cmp
        .descriptor_deltas
        .iter()
        .map(|delta| {
            let dd = PyDict::new(py);
            dd.set_item("left_index", delta.left_index)?;
            dd.set_item("right_index", delta.right_index)?;
            dd.set_item("mw", delta.molecular_weight)?;
            dd.set_item("logp", delta.logp)?;
            dd.set_item("tpsa", delta.tpsa)?;
            dd.set_item("hbd", delta.hbd)?;
            dd.set_item("hba", delta.hba)?;
            dd.set_item("rotatable_bonds", delta.rotatable_bonds)?;
            dd.set_item("qed", delta.qed)?;
            dd.set_item("sa_score", delta.sa_score)?;
            Ok(dd)
        })
        .collect::<PyResult<Vec<_>>>()?;
    d.set_item("descriptor_deltas", deltas)?;
    match &cmp.mcs_smiles {
        Some(s) => d.set_item("mcs_smiles", s)?,
        None => d.set_item("mcs_smiles", py.None())?,
    }
    Ok(d)
}

/// Print a version and accuracy summary — useful for debugging and reporting.
///
/// ```python
/// chematic.doctor()
/// # chematic v1.0.5
/// # Python 3.13.x  |  darwin arm64
/// # ...
/// ```
#[pyfunction]
fn doctor(py: Python<'_>) {
    let ver = env!("CARGO_PKG_VERSION");
    let vi = py.version_info();
    let py_ver = format!("{}.{}.{}", vi.major, vi.minor, vi.patch);

    let platform = std::env::consts::OS;
    let arch = std::env::consts::ARCH;

    println!("chematic v{ver}");
    println!("Python {py_ver}  |  {platform} {arch}");
    println!();
    println!("Descriptor validation (2026-08-23, v0.18.0 vs RDKit 2026.03.4):");
    println!("  MW                     99.82% (4,999-mol, ±0.01 Da)");
    println!("  HBA / HBD / ARC        100.0% (4,999-mol ChEMBL subset)");
    println!("  TPSA                  100.0% (4,999-mol, ±0.1 Å²)");
    println!("  LogP (Crippen)        100.0% (4,999-mol, max Δ = 1.1×10⁻¹³)");
    println!("  Stereocenters         99.96% legacy / 98.6% new CIP (4,999-mol)");
    println!("  CIP R/S/E/Z           99.74% (opt-in accurate engine vs rdCIPLabeler)");
    println!("  ECFP4 throughput      54.7 vs 94.3 µs/mol (v0.18.0/RDKit, 5,000-mol)");
    println!("  WASM bundle           3.30 MB raw / 1.21 MB gzip (v1.0.2 candidate)");
    println!();
    println!("Feature stability:");
    println!("  Stable      SMILES · MW/HBA/HBD/TPSA/LogP · ECFP4/MACCS · SDF/MOL · SMARTS");
    println!("  Stable      Tanimoto · PAINS/Brenk · 2D SVG · QED");
    println!("  Experimental   3D conformer (not RDKit ETKDGv3 equivalent)");
    println!("  Rule-based     pKa · ADMET (screening use only, not clinical)");
    println!("  Partial        IUPAC name generation · pure-Rust InChI (approx.)");
    println!();
    println!("Docs:       https://kent-tokyo.github.io/chematic/");
    println!("Validation: https://github.com/kent-tokyo/chematic/tree/main/validation/");
    println!("Benchmarks: https://github.com/kent-tokyo/chematic/tree/main/benchmarks/");
}

/// Generate a self-contained HTML report for a list of molecules.
///
/// Returns an HTML string. If ``output`` is given, also writes the file.
/// Cards are sorted by QED descending (most drug-like first).
///
/// ```python
/// mols = [chematic.from_smiles(s) for s in ["CC(=O)Oc1ccccc1C(=O)O", "Cn1cnc2c1c(=O)n(C)c(=O)n2C"]]
/// html = chematic.report(mols, names=["aspirin", "caffeine"], output="report.html")
/// ```
#[pyfunction]
#[pyo3(signature = (mols, names=None, title="chematic report", output=None))]
fn report(
    mols: Vec<Mol>,
    names: Option<Vec<Option<String>>>,
    title: &str,
    output: Option<&str>,
) -> PyResult<Report> {
    use chematic_chem::{
        brenk_passes, hbd_count, logp_and_mr, molecular_weight, pains_passes, qed_with_bundle,
        ring_bundle, tpsa,
    };

    // Build (qed, card_html) pairs so we can sort by QED
    let mut cards: Vec<(f64, String)> = mols
        .iter()
        .enumerate()
        .map(|(i, mol)| {
            let m = mol.inner.as_ref();
            let mw = molecular_weight(m);
            let (logp, _) = logp_and_mr(m);
            let tpsa_val = tpsa(m);
            let hbd = hbd_count(m);
            let rb = ring_bundle(m);
            let qed = qed_with_bundle(m, &rb);
            let lip = mw <= 500.0 && hbd <= 5 && rb.hba_count <= 10 && logp <= 5.0;
            let pains_ok = pains_passes(m);
            let brenk_ok = brenk_passes(m);

            let label = names
                .as_ref()
                .and_then(|ns| ns.get(i))
                .and_then(|n| n.as_deref())
                .unwrap_or("");

            let svg = chematic_depict::depict_svg(m);

            let lip_badge = if lip {
                r#"<span class="badge pass">Lipinski ✓</span>"#
            } else {
                r#"<span class="badge fail">Lipinski ✗</span>"#
            };
            let pains_badge = if pains_ok {
                r#"<span class="badge pass">PAINS ✓</span>"#
            } else {
                r#"<span class="badge fail">PAINS ✗</span>"#
            };
            let brenk_badge = if brenk_ok {
                r#"<span class="badge pass">Brenk ✓</span>"#
            } else {
                r#"<span class="badge warn">Brenk ⚠</span>"#
            };

            let card = format!(
                r#"<div class="card"><div class="svg">{svg}</div><div class="name">{label}</div><div class="desc">MW: {mw:.1} Da &nbsp;|&nbsp; LogP: {logp:.2} &nbsp;|&nbsp; TPSA: {tpsa_val:.1} Å²<br>HBD: {hbd} &nbsp;|&nbsp; HBA: {hba} &nbsp;|&nbsp; QED: {qed:.2}</div><div class="badges">{lip_badge}{pains_badge}{brenk_badge}</div></div>"#,
                hba = rb.hba_count,
            );
            (qed, card)
        })
        .collect();

    cards.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));
    let cards_html: String = cards.into_iter().map(|(_, c)| c).collect();
    let n = mols.len();
    let ver = env!("CARGO_PKG_VERSION");

    let html = format!(
        r#"<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{title}</title>
<style>
body{{font-family:system-ui,sans-serif;background:#f8f9fa;padding:24px;margin:0}}
h1{{font-size:1.4rem;color:#333;margin-bottom:4px}}
.meta{{font-size:.85rem;color:#666;margin-bottom:20px}}
.grid{{display:flex;flex-wrap:wrap;gap:16px}}
.card{{background:#fff;border:1px solid #dee2e6;border-radius:8px;padding:12px;width:220px;box-shadow:0 1px 3px rgba(0,0,0,.06)}}
.svg{{width:100%;height:160px;overflow:hidden}}
.svg svg{{width:100%;height:100%}}
.name{{font-weight:600;font-size:.9rem;margin:6px 0 4px;color:#333;min-height:1.1em}}
.desc{{font-size:.78rem;color:#555;line-height:1.8;margin:4px 0}}
.badges{{display:flex;flex-wrap:wrap;gap:4px;margin-top:6px}}
.badge{{font-size:.7rem;padding:2px 7px;border-radius:10px;font-weight:500}}
.pass{{background:#d1e7dd;color:#0a3622}}
.fail{{background:#f8d7da;color:#58151c}}
.warn{{background:#fff3cd;color:#664d03}}
</style>
</head>
<body>
<h1>{title}</h1>
<p class="meta">{n} molecule{plural} &middot; generated by chematic v{ver}</p>
<div class="grid">{cards_html}</div>
</body>
</html>"#,
        plural = if n == 1 { "" } else { "s" },
    );

    let rep = Report { html };
    if let Some(path) = output {
        rep.save(path)?;
    }
    Ok(rep)
}

/// Compare two molecules side-by-side and return a self-contained HTML report.
///
/// Returns a ``Report`` object. In Jupyter, writing ``report`` renders it inline.
///
/// ```python
/// aspirin   = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")
/// ibuprofen = chematic.from_smiles("CC(C)Cc1ccc(CC(C)C(=O)O)cc1")
/// report = chematic.compare(aspirin, ibuprofen, names=("Aspirin", "Ibuprofen"))
/// report.save("compare.html")
/// ```
#[pyfunction]
#[pyo3(signature = (mol1, mol2, names=None, title=None))]
fn compare(mol1: &Mol, mol2: &Mol, names: Option<(String, String)>, title: Option<&str>) -> Report {
    use chematic_chem::{
        brenk_passes, hbd_count, logp_and_mr, molecular_weight, pains_passes, qed_with_bundle,
        ring_bundle, tpsa,
    };

    let m1 = mol1.inner.as_ref();
    let m2 = mol2.inner.as_ref();

    let (name1, name2) = names.unwrap_or_else(|| ("Molecule A".into(), "Molecule B".into()));

    let heading = title
        .map(|t| t.to_string())
        .unwrap_or_else(|| format!("{name1} vs {name2}"));

    let svg1 = chematic_depict::depict_svg(m1);
    let svg2 = chematic_depict::depict_svg(m2);

    let mw1 = molecular_weight(m1);
    let mw2 = molecular_weight(m2);
    let (logp1, _) = logp_and_mr(m1);
    let (logp2, _) = logp_and_mr(m2);
    let tpsa1 = tpsa(m1);
    let tpsa2 = tpsa(m2);
    let hbd1 = hbd_count(m1);
    let hbd2 = hbd_count(m2);
    let rb1 = ring_bundle(m1);
    let rb2 = ring_bundle(m2);
    let qed1 = qed_with_bundle(m1, &rb1);
    let qed2 = qed_with_bundle(m2, &rb2);
    let lip1 = mw1 <= 500.0 && hbd1 <= 5 && rb1.hba_count <= 10 && logp1 <= 5.0;
    let lip2 = mw2 <= 500.0 && hbd2 <= 5 && rb2.hba_count <= 10 && logp2 <= 5.0;
    let pains1 = pains_passes(m1);
    let pains2 = pains_passes(m2);
    let brenk1 = brenk_passes(m1);
    let brenk2 = brenk_passes(m2);

    // MCS common atoms (reuse diff logic)
    let config = chematic_smarts::McsConfig::default();
    let qmol = chematic_smarts::find_mcs_with_config(&[m1, m2], &config);
    let common = qmol.atoms.len();

    // Delta summary (same logic as Mol::diff)
    let elem_parts: Vec<String> = {
        use std::collections::BTreeMap;
        let mut c1: BTreeMap<&'static str, i32> = BTreeMap::new();
        let mut c2: BTreeMap<&'static str, i32> = BTreeMap::new();
        for i in 0..m1.atom_count() {
            *c1.entry(m1.atom(chematic_core::AtomIdx(i as u32)).element.symbol())
                .or_insert(0) += 1;
        }
        for i in 0..m2.atom_count() {
            *c2.entry(m2.atom(chematic_core::AtomIdx(i as u32)).element.symbol())
                .or_insert(0) += 1;
        }
        let all: std::collections::BTreeSet<_> = c1.keys().chain(c2.keys()).copied().collect();
        all.iter()
            .filter_map(|e| {
                let d = c2.get(e).copied().unwrap_or(0) - c1.get(e).copied().unwrap_or(0);
                if d != 0 {
                    Some(if d > 0 {
                        format!("+{d}{e}")
                    } else {
                        format!("{d}{e}")
                    })
                } else {
                    None
                }
            })
            .collect()
    };
    let elem_str = if elem_parts.is_empty() {
        "Same elemental composition".into()
    } else {
        elem_parts.join(", ")
    };
    let summary = format!(
        "{}. \u{0394}LogP {:+.2}, \u{0394}TPSA {:+.1} \u{00c5}\u{00b2}, \u{0394}MW {:+.1} Da.",
        elem_str,
        logp2 - logp1,
        tpsa2 - tpsa1,
        mw2 - mw1,
    );

    fn delta_class(d: f64) -> &'static str {
        if d > 0.0 {
            "pos"
        } else if d < 0.0 {
            "neg"
        } else {
            ""
        }
    }
    fn flag(v: bool, ok: &str, fail: &str) -> String {
        if v {
            format!(r#"<span class="pass">{ok}</span>"#)
        } else {
            format!(r#"<span class="fail">{fail}</span>"#)
        }
    }
    fn warn_flag(v: bool) -> String {
        if v {
            r#"<span class="pass">✓</span>"#.into()
        } else {
            r#"<span class="warn">⚠</span>"#.into()
        }
    }

    let ver = env!("CARGO_PKG_VERSION");
    let html = format!(
        r#"<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{heading}</title>
<style>
body{{font-family:system-ui,sans-serif;background:#f8f9fa;padding:24px;margin:0}}
h1{{font-size:1.4rem;color:#333;margin-bottom:4px}}
.meta{{font-size:.85rem;color:#666;margin-bottom:6px}}
.summary{{font-size:.85rem;color:#444;background:#fff;border:1px solid #dee2e6;border-radius:6px;padding:8px 12px;margin-bottom:20px;display:inline-block}}
table{{border-collapse:collapse;background:#fff;border-radius:8px;overflow:hidden;box-shadow:0 1px 3px rgba(0,0,0,.06)}}
th,td{{padding:10px 16px;text-align:left;border-bottom:1px solid #f0f0f0;font-size:.88rem}}
th{{background:#f8f9fa;font-weight:600;color:#555}}
td.num{{text-align:right;font-variant-numeric:tabular-nums}}
td.delta{{text-align:right;font-size:.8rem;font-weight:600}}
.pos{{color:#0a6640}}
.neg{{color:#8b1c1c}}
.pass{{color:#0a6640;font-weight:600}}
.fail{{color:#8b1c1c;font-weight:600}}
.warn{{color:#7d5a00;font-weight:600}}
.svg-cell svg{{width:180px;height:140px}}
</style>
</head>
<body>
<h1>{heading}</h1>
<p class="meta">Common scaffold: {common} atoms &middot; generated by chematic v{ver}</p>
<p class="summary">{summary}</p>
<table>
<tr><th>Property</th><th>{name1}</th><th>{name2}</th><th>Delta</th></tr>
<tr><td>Structure</td>
    <td class="svg-cell">{svg1}</td>
    <td class="svg-cell">{svg2}</td>
    <td></td></tr>
<tr><td>MW (Da)</td>
    <td class="num">{mw1:.1}</td><td class="num">{mw2:.1}</td>
    <td class="delta {dc_mw}">{dmw:+.1}</td></tr>
<tr><td>LogP</td>
    <td class="num">{logp1:.2}</td><td class="num">{logp2:.2}</td>
    <td class="delta {dc_lp}">{dlp:+.2}</td></tr>
<tr><td>TPSA (Å²)</td>
    <td class="num">{tpsa1:.1}</td><td class="num">{tpsa2:.1}</td>
    <td class="delta {dc_tp}">{dtp:+.1}</td></tr>
<tr><td>HBD</td>
    <td class="num">{hbd1}</td><td class="num">{hbd2}</td>
    <td class="delta {dc_hbd}">{dhbd:+}</td></tr>
<tr><td>HBA</td>
    <td class="num">{hba1}</td><td class="num">{hba2}</td>
    <td class="delta {dc_hba}">{dhba:+}</td></tr>
<tr><td>QED</td>
    <td class="num">{qed1:.2}</td><td class="num">{qed2:.2}</td>
    <td class="delta {dc_qed}">{dqed:+.2}</td></tr>
<tr><td>Lipinski</td>
    <td>{lip1_s}</td><td>{lip2_s}</td><td></td></tr>
<tr><td>PAINS</td>
    <td>{pains1_s}</td><td>{pains2_s}</td><td></td></tr>
<tr><td>Brenk</td>
    <td>{brenk1_s}</td><td>{brenk2_s}</td><td></td></tr>
</table>
</body>
</html>"#,
        dc_mw = delta_class(mw2 - mw1),
        dmw = mw2 - mw1,
        dc_lp = delta_class(logp2 - logp1),
        dlp = logp2 - logp1,
        dc_tp = delta_class(tpsa2 - tpsa1),
        dtp = tpsa2 - tpsa1,
        dc_hbd = delta_class((hbd2 as f64) - (hbd1 as f64)),
        dhbd = hbd2 as i32 - hbd1 as i32,
        dc_hba = delta_class((rb2.hba_count as f64) - (rb1.hba_count as f64)),
        dhba = rb2.hba_count as i32 - rb1.hba_count as i32,
        hba1 = rb1.hba_count,
        hba2 = rb2.hba_count,
        dc_qed = delta_class(qed2 - qed1),
        dqed = qed2 - qed1,
        lip1_s = flag(lip1, "✓ Lipinski", "✗ Lipinski"),
        lip2_s = flag(lip2, "✓ Lipinski", "✗ Lipinski"),
        pains1_s = flag(pains1, "✓ PAINS", "✗ PAINS"),
        pains2_s = flag(pains2, "✓ PAINS", "✗ PAINS"),
        brenk1_s = warn_flag(brenk1),
        brenk2_s = warn_flag(brenk2),
    );

    Report { html }
}

// ---------------------------------------------------------------------------
// Register
// ---------------------------------------------------------------------------

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Report>()?;
    m.add_class::<Report>()?;
    m.add_function(wrap_pyfunction!(molecule_report, m)?)?;
    m.add_function(wrap_pyfunction!(screen_smiles, m)?)?;
    m.add_function(wrap_pyfunction!(compare_molecules, m)?)?;
    m.add_function(wrap_pyfunction!(doctor, m)?)?;
    m.add_function(wrap_pyfunction!(report, m)?)?;
    m.add_function(wrap_pyfunction!(compare, m)?)?;
    Ok(())
}
