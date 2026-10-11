use ndarray::{Array1, Array2};
use numpy::{IntoPyArray, PyArray1, PyArray2};
use pyo3::prelude::*;
use pyo3::types::PyDict;
use rayon::prelude::*;
use std::sync::Arc;

use crate::Mol;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn parse_parallel(smiles: &[String]) -> Vec<Option<chematic_core::Molecule>> {
    smiles
        .par_iter()
        .map(|s| chematic_smiles::parse(s).ok())
        .collect()
}

fn bitvec2048_to_bits(fp: &chematic_fp::bitvec::BitVec2048) -> Vec<u8> {
    (0..2048usize).map(|i| u8::from(fp.get(i))).collect()
}

/// Internal row representation used while computing `bulk.descriptors_array`.
/// Keeping this at module scope lets the expensive row computation and the
/// NumPy materialization live in separate, testable units.
#[derive(Default)]
struct DescriptorRow {
    // floats (int fields cast to f64 for uniform storage)
    mw: f64,
    exact_mass: f64,
    tpsa: f64,
    logp: f64,
    mr: f64,
    hbd: f64,
    hba: f64,
    rb: f64,
    hac: f64,
    rc: f64,
    arc: f64,
    nh: f64,
    nsc: f64,
    nsp: f64,
    nbh: f64,
    fsp3: f64,
    qed: f64,
    sa: f64,
    fc: f64,
    asa: f64,
    bertz: f64,
    wi: f64,
    k1: f64,
    k2: f64,
    k3: f64,
    c0: f64,
    c1: f64,
    c2: f64,
    c3: f64,
    c4: f64,
    c0v: f64,
    c1v: f64,
    c2v: f64,
    c3v: f64,
    c4v: f64,
    n_ah: f64,
    n_alh: f64,
    n_sr: f64,
    n_ar: f64,
    n_usc: f64,
    sum_e: f64,
    max_e: f64,
    min_e: f64,
    bbb: f64,
    caco: f64,
    herg: f64,
    cyp: f64,
    schultz: f64,
    gutman: f64,
    vabc: f64,
    grav: f64,
    // booleans
    lip: bool,
    veb: bool,
    egan: bool,
    ghose: bool,
    reos: bool,
    pains: bool,
    bbp: bool,
    // optional floats
    pka_acid: Option<f64>,
    pka_base: Option<f64>,
}

const VALID_DESCRIPTOR_COLUMNS: &[&str] = &[
    "mw",
    "exact_mass",
    "tpsa",
    "logp",
    "molar_refractivity",
    "hbd",
    "hba",
    "rotatable_bonds",
    "heavy_atoms",
    "ring_count",
    "aromatic_ring_count",
    "num_heteroatoms",
    "num_stereocenters",
    "num_spiro_atoms",
    "num_bridgehead_atoms",
    "fsp3",
    "qed",
    "sa_score",
    "formal_charge",
    "labute_asa",
    "bertz_ct",
    "wiener_index",
    "kappa1",
    "kappa2",
    "kappa3",
    "chi0",
    "chi1",
    "chi2",
    "chi3",
    "chi4",
    "chi0v",
    "chi1v",
    "chi2v",
    "chi3v",
    "chi4v",
    "num_aromatic_heterocycles",
    "num_aliphatic_heterocycles",
    "num_saturated_rings",
    "num_aliphatic_rings",
    "num_unspecified_stereocenters",
    "sum_estate",
    "max_estate",
    "min_estate",
    "lipinski_passes",
    "veber_passes",
    "egan_passes",
    "ghose_passes",
    "reos_passes",
    "pains_passes",
    "bbb_passes",
    "bbb_score",
    "caco2",
    "herg_risk",
    "cyp3a4_risk",
    "pka_acid",
    "pka_base",
    "schultz_mti",
    "gutman_mti",
    "vabc",
    "gravitational_index",
];

struct DescriptorColumnPlan<'a> {
    columns: &'a [String],
    need_qed: bool,
    need_sa: bool,
    need_pka: bool,
    need_logp: bool,
    need_mr: bool,
    need_tpsa: bool,
    need_hbd: bool,
    need_ring_bundle: bool,
    need_kappa: bool,
    need_chi: bool,
    need_estate: bool,
    need_mw: bool,
    need_hac: bool,
    need_fc: bool,
}

impl<'a> DescriptorColumnPlan<'a> {
    fn new(columns: &'a [String]) -> PyResult<Self> {
        for column in columns {
            if !VALID_DESCRIPTOR_COLUMNS.contains(&column.as_str()) {
                return Err(pyo3::exceptions::PyValueError::new_err(format!(
                    "unknown column: {column:?}"
                )));
            }
        }
        let wants_any = |names: &[&str]| names.iter().any(|name| wants(columns, name));
        Ok(Self {
            columns,
            need_qed: wants(columns, "qed"),
            need_sa: wants(columns, "sa_score"),
            need_pka: wants_any(&["pka_acid", "pka_base"]),
            need_logp: wants_any(&[
                "logp",
                "molar_refractivity",
                "qed",
                "bbb_score",
                "caco2",
                "herg_risk",
                "cyp3a4_risk",
                "lipinski_passes",
                "egan_passes",
                "ghose_passes",
                "reos_passes",
            ]),
            need_mr: wants_any(&["molar_refractivity", "ghose_passes"]),
            need_tpsa: wants_any(&[
                "tpsa",
                "qed",
                "bbb_score",
                "bbb_passes",
                "caco2",
                "veber_passes",
                "egan_passes",
            ]),
            need_hbd: wants_any(&["hbd", "qed", "bbb_passes", "lipinski_passes", "reos_passes"]),
            need_ring_bundle: wants_any(&[
                "hba",
                "rotatable_bonds",
                "ring_count",
                "aromatic_ring_count",
                "num_spiro_atoms",
                "num_bridgehead_atoms",
                "num_aromatic_heterocycles",
                "num_aliphatic_heterocycles",
                "num_saturated_rings",
                "num_aliphatic_rings",
                "qed",
                "sa_score",
                "cyp3a4_risk",
                "lipinski_passes",
                "reos_passes",
            ]),
            need_kappa: wants_any(&["kappa1", "kappa2", "kappa3"]),
            need_chi: wants_any(&[
                "chi0", "chi1", "chi2", "chi3", "chi4", "chi0v", "chi1v", "chi2v", "chi3v", "chi4v",
            ]),
            need_estate: wants_any(&["sum_estate", "max_estate", "min_estate"]),
            need_mw: wants_any(&[
                "mw",
                "bbb_passes",
                "bbb_score",
                "herg_risk",
                "cyp3a4_risk",
                "lipinski_passes",
                "ghose_passes",
                "reos_passes",
            ]),
            need_hac: wants_any(&["heavy_atoms", "ghose_passes", "reos_passes"]),
            need_fc: wants_any(&["formal_charge", "reos_passes"]),
        })
    }

    fn wants(&self, name: &str) -> bool {
        wants(self.columns, name)
    }
}

fn wants(columns: &[String], name: &str) -> bool {
    columns.iter().any(|column| column == name)
}

struct SelectedDescriptorContext<'a, 'p> {
    molecule: &'a chematic_core::Molecule,
    plan: &'p DescriptorColumnPlan<'p>,
    rings: Option<chematic_chem::RingBundle>,
    topology: Option<chematic_chem::TopologyBundle>,
    pka: (Option<f64>, Option<f64>),
    logp_mr: (f64, f64),
    kappa: (f64, f64, f64),
    chi: (f64, f64, f64, f64, f64, f64, f64, f64, f64, f64),
    estate: (f64, f64, f64),
    mw: f64,
    tpsa: f64,
    hbd: usize,
    hba: usize,
    rotatable_bonds: usize,
    heavy_atoms: usize,
    formal_charge: i32,
}

impl<'a, 'p> SelectedDescriptorContext<'a, 'p> {
    fn new(molecule: &'a chematic_core::Molecule, plan: &'p DescriptorColumnPlan<'p>) -> Self {
        let rings = plan
            .need_ring_bundle
            .then(|| chematic_chem::ring_bundle(molecule));
        let topology_requests = [plan.need_kappa, plan.need_chi, plan.wants("wiener_index")]
            .into_iter()
            .filter(|requested| *requested)
            .count();
        let topology = (topology_requests >= 2).then(|| chematic_chem::topology_bundle(molecule));
        let pka = if plan.need_pka || plan.wants("herg_risk") {
            chematic_chem::pka_both(molecule)
        } else {
            (None, None)
        };
        let logp_mr = if plan.need_logp || plan.need_mr {
            chematic_chem::logp_and_mr(molecule)
        } else {
            (0.0, 0.0)
        };
        let kappa = if plan.need_kappa {
            topology
                .as_ref()
                .map(|value| value.kappa)
                .unwrap_or_else(|| chematic_chem::kappa_all(molecule))
        } else {
            (0.0, 0.0, 0.0)
        };
        let chi = if plan.need_chi {
            topology
                .as_ref()
                .map(|value| value.chi)
                .unwrap_or_else(|| chematic_chem::chi_all(molecule))
        } else {
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        };
        let estate = if plan.need_estate {
            chematic_chem::estate_all(molecule)
        } else {
            (0.0, 0.0, 0.0)
        };
        let mw = if plan.need_mw {
            chematic_chem::molecular_weight(molecule)
        } else {
            0.0
        };
        let tpsa = if plan.need_tpsa {
            chematic_chem::tpsa(molecule)
        } else {
            0.0
        };
        let hbd = if plan.need_hbd {
            chematic_chem::hbd_count(molecule)
        } else {
            0
        };
        let hba = rings.as_ref().map_or(0, |value| value.hba_count);
        let rotatable_bonds = rings.as_ref().map_or(0, |value| value.rotatable_bond_count);
        let heavy_atoms = if plan.need_hac {
            chematic_chem::heavy_atom_count(molecule)
        } else {
            0
        };
        let formal_charge = if plan.need_fc {
            chematic_chem::formal_charge_sum(molecule)
        } else {
            0
        };
        Self {
            molecule,
            plan,
            rings,
            topology,
            pka,
            logp_mr,
            kappa,
            chi,
            estate,
            mw,
            tpsa,
            hbd,
            hba,
            rotatable_bonds,
            heavy_atoms,
            formal_charge,
        }
    }

    fn into_row(self) -> DescriptorRow {
        let mut row = DescriptorRow::default();
        self.fill_basic(&mut row);
        self.fill_topology(&mut row);
        self.fill_admet(&mut row);
        self.fill_filters(&mut row);
        row
    }

    fn fill_basic(&self, row: &mut DescriptorRow) {
        row.mw = self.mw;
        row.exact_mass = if self.plan.wants("exact_mass") {
            chematic_chem::exact_mass(self.molecule)
        } else {
            0.0
        };
        row.tpsa = self.tpsa;
        row.logp = self.logp_mr.0;
        row.mr = self.logp_mr.1;
        row.hbd = self.hbd as f64;
        row.hba = self.hba as f64;
        row.rb = self.rotatable_bonds as f64;
        row.hac = self.heavy_atoms as f64;
        row.fsp3 = if self.plan.wants("fsp3") {
            chematic_chem::fsp3(self.molecule)
        } else {
            0.0
        };
        row.qed = if self.plan.need_qed {
            chematic_chem::qed_with_bundle(
                self.molecule,
                self.rings.as_ref().expect("QED requires ring bundle"),
            )
        } else {
            0.0
        };
        row.sa = if self.plan.need_sa {
            chematic_chem::sa_score_with_bundle(
                self.molecule,
                self.rings.as_ref().expect("SA requires ring bundle"),
            )
        } else {
            0.0
        };
        row.fc = self.formal_charge as f64;
        row.asa = if self.plan.wants("labute_asa") {
            chematic_chem::labute_asa(self.molecule)
        } else {
            0.0
        };
        row.bertz = if self.plan.wants("bertz_ct") {
            chematic_chem::bertz_ct(self.molecule)
        } else {
            0.0
        };
        row.pka_acid = self.pka.0;
        row.pka_base = self.pka.1;
    }

    fn fill_topology(&self, row: &mut DescriptorRow) {
        let rings = self.rings.as_ref();
        row.rc = rings.map_or(0, |value| value.ring_count) as f64;
        row.arc = rings.map_or(0, |value| value.aromatic_ring_count) as f64;
        row.nh = if self.plan.wants("num_heteroatoms") {
            chematic_chem::num_heteroatoms(self.molecule)
        } else {
            0
        } as f64;
        row.nsc = if self.plan.wants("num_stereocenters") {
            chematic_chem::num_stereocenters(self.molecule)
        } else {
            0
        } as f64;
        row.nsp = rings.map_or(0, |value| value.num_spiro_atoms) as f64;
        row.nbh = rings.map_or(0, |value| value.num_bridgehead_atoms) as f64;
        (row.k1, row.k2, row.k3) = self.kappa;
        (
            row.c0, row.c1, row.c2, row.c3, row.c4, row.c0v, row.c1v, row.c2v, row.c3v, row.c4v,
        ) = self.chi;
        row.wi = if self.plan.wants("wiener_index") {
            self.topology
                .as_ref()
                .map(|value| value.wiener)
                .unwrap_or_else(|| chematic_chem::wiener_index(self.molecule))
        } else {
            0.0
        };
        row.n_ah = rings.map_or(0, |value| value.num_aromatic_heterocycles) as f64;
        row.n_alh = rings.map_or(0, |value| value.num_aliphatic_heterocycles) as f64;
        row.n_sr = rings.map_or(0, |value| value.num_saturated_rings) as f64;
        row.n_ar = rings.map_or(0, |value| value.num_aliphatic_rings) as f64;
        row.n_usc = if self.plan.wants("num_unspecified_stereocenters") {
            chematic_chem::num_unspecified_stereocenters(self.molecule)
        } else {
            0
        } as f64;
        (row.sum_e, row.max_e, row.min_e) = self.estate;
    }

    fn fill_admet(&self, row: &mut DescriptorRow) {
        row.bbb = if self.plan.wants("bbb_score") {
            chematic_chem::bbb_score_from_parts(self.tpsa, self.logp_mr.0)
        } else {
            0.0
        };
        row.caco = if self.plan.wants("caco2") {
            chematic_chem::caco2_precomputed(self.tpsa, self.logp_mr.0)
        } else {
            0.0
        };
        row.herg = if self.plan.wants("herg_risk") {
            chematic_chem::herg_risk_precomputed(self.molecule, self.logp_mr.0, self.mw)
        } else {
            0.0
        };
        row.cyp = if self.plan.wants("cyp3a4_risk") {
            {
                chematic_chem::cyp3a4_precomputed(
                    self.mw,
                    self.logp_mr.0,
                    chematic_chem::num_aromatic_heterocycles(self.molecule),
                    self.hba,
                )
            }
        } else {
            0.0
        };
        row.schultz = if self.plan.wants("schultz_mti") {
            chematic_chem::schultz_mti(self.molecule) as f64
        } else {
            0.0
        };
        row.gutman = if self.plan.wants("gutman_mti") {
            chematic_chem::gutman_mti(self.molecule) as f64
        } else {
            0.0
        };
        row.vabc = if self.plan.wants("vabc") {
            chematic_chem::vabc(self.molecule)
        } else {
            0.0
        };
        row.grav = if self.plan.wants("gravitational_index") {
            chematic_chem::gravitational_index(self.molecule)
        } else {
            0.0
        };
    }

    fn fill_filters(&self, row: &mut DescriptorRow) {
        row.lip = self.plan.wants("lipinski_passes")
            && self.mw <= 500.0
            && self.hbd <= 5
            && self.hba <= 10
            && self.logp_mr.0 <= 5.0;
        row.veb =
            self.plan.wants("veber_passes") && self.tpsa <= 140.0 && self.rotatable_bonds <= 10;
        row.egan = self.plan.wants("egan_passes") && self.tpsa <= 131.6 && self.logp_mr.0 <= 5.88;
        row.ghose = self.plan.wants("ghose_passes")
            && (160.0..=480.0).contains(&self.mw)
            && (-0.4..=5.6).contains(&self.logp_mr.0)
            && (20.0..=70.0).contains(&(self.heavy_atoms as f64))
            && (40.0..=130.0).contains(&self.logp_mr.1);
        row.reos = self.plan.wants("reos_passes")
            && (200.0..=500.0).contains(&self.mw)
            && (-5.0..=5.0).contains(&self.logp_mr.0)
            && (0..=5).contains(&self.hbd)
            && (0..=10).contains(&self.hba)
            && (-2..=2).contains(&self.formal_charge)
            && (0..=8).contains(&self.rotatable_bonds)
            && (15..=50).contains(&self.heavy_atoms);
        row.pains = self.plan.wants("pains_passes") && chematic_chem::pains_passes(self.molecule);
        row.bbp =
            self.plan.wants("bbb_passes") && self.tpsa < 90.0 && self.mw < 400.0 && self.hbd <= 3;
    }
}

struct DescriptorRecord {
    row: DescriptorRow,
    hbd: usize,
    hba: usize,
    rotatable_bonds: usize,
    heavy_atoms: usize,
    ring_count: usize,
    aromatic_ring_count: usize,
    num_heteroatoms: usize,
    num_stereocenters: usize,
    num_spiro_atoms: usize,
    num_bridgehead_atoms: usize,
    num_aromatic_heterocycles: usize,
    num_aliphatic_heterocycles: usize,
    num_saturated_rings: usize,
    num_aliphatic_rings: usize,
    num_unspecified_stereocenters: usize,
    formal_charge: i32,
    schultz_mti: u64,
    gutman_mti: u64,
}

fn compute_descriptor_record(molecule: &chematic_core::Molecule) -> DescriptorRecord {
    let rings = chematic_chem::ring_bundle(molecule);
    let topology = chematic_chem::topology_bundle(molecule);
    let (pka_acid, pka_base) = chematic_chem::pka_both(molecule);
    let (logp, mr) = chematic_chem::logp_and_mr(molecule);
    let (k1, k2, k3) = topology.kappa;
    let (c0, c1, c2, c3, c4, c0v, c1v, c2v, c3v, c4v) = topology.chi;
    let (sum_e, max_e, min_e) = chematic_chem::estate_all(molecule);
    let hbd = chematic_chem::hbd_count(molecule);
    let heavy_atoms = chematic_chem::heavy_atom_count(molecule);
    let formal_charge = chematic_chem::formal_charge_sum(molecule);
    DescriptorRecord {
        row: DescriptorRow {
            mw: chematic_chem::molecular_weight(molecule),
            exact_mass: chematic_chem::exact_mass(molecule),
            tpsa: chematic_chem::tpsa(molecule),
            logp,
            mr,
            hbd: hbd as f64,
            hba: rings.hba_count as f64,
            rb: rings.rotatable_bond_count as f64,
            hac: heavy_atoms as f64,
            rc: rings.ring_count as f64,
            arc: rings.aromatic_ring_count as f64,
            nh: chematic_chem::num_heteroatoms(molecule) as f64,
            nsc: chematic_chem::num_stereocenters(molecule) as f64,
            nsp: rings.num_spiro_atoms as f64,
            nbh: rings.num_bridgehead_atoms as f64,
            fsp3: chematic_chem::fsp3(molecule),
            qed: chematic_chem::qed(molecule),
            sa: chematic_chem::sa_score(molecule),
            fc: formal_charge as f64,
            asa: chematic_chem::labute_asa(molecule),
            bertz: chematic_chem::bertz_ct(molecule),
            wi: chematic_chem::wiener_index(molecule),
            k1,
            k2,
            k3,
            c0,
            c1,
            c2,
            c3,
            c4,
            c0v,
            c1v,
            c2v,
            c3v,
            c4v,
            n_ah: rings.num_aromatic_heterocycles as f64,
            n_alh: rings.num_aliphatic_heterocycles as f64,
            n_sr: rings.num_saturated_rings as f64,
            n_ar: rings.num_aliphatic_rings as f64,
            n_usc: chematic_chem::num_unspecified_stereocenters(molecule) as f64,
            sum_e,
            max_e,
            min_e,
            bbb: chematic_chem::bbb_score(molecule),
            caco: chematic_chem::caco2_permeability(molecule),
            herg: chematic_chem::herg_risk_score(molecule),
            cyp: chematic_chem::cyp3a4_inhibition_risk(molecule),
            schultz: chematic_chem::schultz_mti(molecule) as f64,
            gutman: chematic_chem::gutman_mti(molecule) as f64,
            vabc: chematic_chem::vabc(molecule),
            grav: chematic_chem::gravitational_index(molecule),
            lip: chematic_chem::lipinski_passes(molecule),
            veb: chematic_chem::veber_passes(molecule),
            egan: chematic_chem::egan_passes(molecule),
            ghose: chematic_chem::ghose_passes(molecule),
            reos: chematic_chem::reos_passes(molecule),
            pains: chematic_chem::pains_passes(molecule),
            bbp: chematic_chem::bbb_passes(molecule),
            pka_acid,
            pka_base,
        },
        hbd,
        hba: rings.hba_count,
        rotatable_bonds: rings.rotatable_bond_count,
        heavy_atoms,
        ring_count: rings.ring_count,
        aromatic_ring_count: rings.aromatic_ring_count,
        num_heteroatoms: chematic_chem::num_heteroatoms(molecule),
        num_stereocenters: chematic_chem::num_stereocenters(molecule),
        num_spiro_atoms: rings.num_spiro_atoms,
        num_bridgehead_atoms: rings.num_bridgehead_atoms,
        num_aromatic_heterocycles: rings.num_aromatic_heterocycles,
        num_aliphatic_heterocycles: rings.num_aliphatic_heterocycles,
        num_saturated_rings: rings.num_saturated_rings,
        num_aliphatic_rings: rings.num_aliphatic_rings,
        num_unspecified_stereocenters: chematic_chem::num_unspecified_stereocenters(molecule),
        formal_charge,
        schultz_mti: chematic_chem::schultz_mti(molecule),
        gutman_mti: chematic_chem::gutman_mti(molecule),
    }
}

fn descriptor_record_to_pydict<'py>(
    py: Python<'py>,
    record: DescriptorRecord,
) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    set_descriptor_scalar_fields(&d, &record)?;
    set_descriptor_topology_fields(&d, &record)?;
    set_descriptor_screening_fields(py, &d, &record)?;
    Ok(d)
}

fn set_descriptor_scalar_fields(d: &Bound<'_, PyDict>, record: &DescriptorRecord) -> PyResult<()> {
    let r = &record.row;
    for (name, value) in [
        ("mw", r.mw),
        ("exact_mass", r.exact_mass),
        ("tpsa", r.tpsa),
        ("logp", r.logp),
        ("molar_refractivity", r.mr),
        ("fsp3", r.fsp3),
        ("qed", r.qed),
        ("sa_score", r.sa),
        ("labute_asa", r.asa),
        ("bertz_ct", r.bertz),
        ("wiener_index", r.wi),
    ] {
        d.set_item(name, value)?;
    }
    d.set_item("hbd", record.hbd)?;
    d.set_item("hba", record.hba)?;
    d.set_item("rotatable_bonds", record.rotatable_bonds)?;
    d.set_item("heavy_atoms", record.heavy_atoms)?;
    d.set_item("formal_charge", record.formal_charge)?;
    Ok(())
}

fn set_descriptor_topology_fields(
    d: &Bound<'_, PyDict>,
    record: &DescriptorRecord,
) -> PyResult<()> {
    let r = &record.row;
    d.set_item("ring_count", record.ring_count)?;
    d.set_item("aromatic_ring_count", record.aromatic_ring_count)?;
    d.set_item("num_heteroatoms", record.num_heteroatoms)?;
    d.set_item("num_stereocenters", record.num_stereocenters)?;
    d.set_item("num_spiro_atoms", record.num_spiro_atoms)?;
    d.set_item("num_bridgehead_atoms", record.num_bridgehead_atoms)?;
    d.set_item(
        "num_aromatic_heterocycles",
        record.num_aromatic_heterocycles,
    )?;
    d.set_item(
        "num_aliphatic_heterocycles",
        record.num_aliphatic_heterocycles,
    )?;
    d.set_item("num_saturated_rings", record.num_saturated_rings)?;
    d.set_item("num_aliphatic_rings", record.num_aliphatic_rings)?;
    d.set_item(
        "num_unspecified_stereocenters",
        record.num_unspecified_stereocenters,
    )?;
    for (name, value) in [
        ("kappa1", r.k1),
        ("kappa2", r.k2),
        ("kappa3", r.k3),
        ("chi0", r.c0),
        ("chi1", r.c1),
        ("chi2", r.c2),
        ("chi3", r.c3),
        ("chi4", r.c4),
        ("chi0v", r.c0v),
        ("chi1v", r.c1v),
        ("chi2v", r.c2v),
        ("chi3v", r.c3v),
        ("chi4v", r.c4v),
        ("sum_estate", r.sum_e),
        ("max_estate", r.max_e),
        ("min_estate", r.min_e),
    ] {
        d.set_item(name, value)?;
    }
    Ok(())
}

fn set_descriptor_screening_fields(
    py: Python<'_>,
    d: &Bound<'_, PyDict>,
    record: &DescriptorRecord,
) -> PyResult<()> {
    let r = &record.row;
    for (name, value) in [
        ("lipinski_passes", r.lip),
        ("veber_passes", r.veb),
        ("egan_passes", r.egan),
        ("ghose_passes", r.ghose),
        ("reos_passes", r.reos),
        ("pains_passes", r.pains),
        ("bbb_passes", r.bbp),
    ] {
        d.set_item(name, value)?;
    }
    for (name, value) in [
        ("bbb_score", r.bbb),
        ("caco2", r.caco),
        ("herg_risk", r.herg),
        ("cyp3a4_risk", r.cyp),
        ("vabc", r.vabc),
        ("gravitational_index", r.grav),
    ] {
        d.set_item(name, value)?;
    }
    match r.pka_acid {
        Some(value) => d.set_item("pka_acid", value)?,
        None => d.set_item("pka_acid", py.None())?,
    }
    match r.pka_base {
        Some(value) => d.set_item("pka_base", value)?,
        None => d.set_item("pka_base", py.None())?,
    }
    d.set_item("schultz_mti", record.schultz_mti)?;
    d.set_item("gutman_mti", record.gutman_mti)?;
    Ok(())
}

fn materialize_descriptor_columns<'py>(
    py: Python<'py>,
    rows: &[DescriptorRow],
    columns: &[String],
) -> PyResult<Bound<'py, PyDict>> {
    let out = PyDict::new(py);

    macro_rules! fcol {
        ($name:literal, $field:ident) => {
            if columns.iter().any(|column| column == $name) {
                let arr = Array1::from(rows.iter().map(|r| r.$field).collect::<Vec<f64>>());
                out.set_item($name, arr.into_pyarray(py))?;
            }
        };
    }
    macro_rules! bcol {
        ($name:literal, $field:ident) => {
            if columns.iter().any(|column| column == $name) {
                let arr = Array1::from(rows.iter().map(|r| r.$field).collect::<Vec<bool>>());
                out.set_item($name, arr.into_pyarray(py))?;
            }
        };
    }
    macro_rules! ocol {
        ($name:literal, $field:ident) => {
            if columns.iter().any(|column| column == $name) {
                let arr = Array1::from(
                    rows.iter()
                        .map(|r| r.$field.unwrap_or(f64::NAN))
                        .collect::<Vec<f64>>(),
                );
                out.set_item($name, arr.into_pyarray(py))?;
            }
        };
    }

    fcol!("mw", mw);
    fcol!("exact_mass", exact_mass);
    fcol!("tpsa", tpsa);
    fcol!("logp", logp);
    fcol!("molar_refractivity", mr);
    fcol!("hbd", hbd);
    fcol!("hba", hba);
    fcol!("rotatable_bonds", rb);
    fcol!("heavy_atoms", hac);
    fcol!("ring_count", rc);
    fcol!("aromatic_ring_count", arc);
    fcol!("num_heteroatoms", nh);
    fcol!("num_stereocenters", nsc);
    fcol!("num_spiro_atoms", nsp);
    fcol!("num_bridgehead_atoms", nbh);
    fcol!("fsp3", fsp3);
    fcol!("qed", qed);
    fcol!("sa_score", sa);
    fcol!("formal_charge", fc);
    fcol!("labute_asa", asa);
    fcol!("bertz_ct", bertz);
    fcol!("wiener_index", wi);
    fcol!("kappa1", k1);
    fcol!("kappa2", k2);
    fcol!("kappa3", k3);
    fcol!("chi0", c0);
    fcol!("chi1", c1);
    fcol!("chi2", c2);
    fcol!("chi3", c3);
    fcol!("chi4", c4);
    fcol!("chi0v", c0v);
    fcol!("chi1v", c1v);
    fcol!("chi2v", c2v);
    fcol!("chi3v", c3v);
    fcol!("chi4v", c4v);
    fcol!("num_aromatic_heterocycles", n_ah);
    fcol!("num_aliphatic_heterocycles", n_alh);
    fcol!("num_saturated_rings", n_sr);
    fcol!("num_aliphatic_rings", n_ar);
    fcol!("num_unspecified_stereocenters", n_usc);
    fcol!("sum_estate", sum_e);
    fcol!("max_estate", max_e);
    fcol!("min_estate", min_e);
    fcol!("bbb_score", bbb);
    fcol!("caco2", caco);
    fcol!("herg_risk", herg);
    fcol!("cyp3a4_risk", cyp);
    fcol!("schultz_mti", schultz);
    fcol!("gutman_mti", gutman);
    fcol!("vabc", vabc);
    fcol!("gravitational_index", grav);
    bcol!("lipinski_passes", lip);
    bcol!("veber_passes", veb);
    bcol!("egan_passes", egan);
    bcol!("ghose_passes", ghose);
    bcol!("reos_passes", reos);
    bcol!("pains_passes", pains);
    bcol!("bbb_passes", bbp);
    ocol!("pka_acid", pka_acid);
    ocol!("pka_base", pka_base);

    Ok(out)
}

// ---------------------------------------------------------------------------
// bulk.parse — parallel SMILES parsing
// ---------------------------------------------------------------------------

/// Parse a list of SMILES strings in parallel.
///
/// Returns a list of ``Mol`` objects (or ``None`` for invalid SMILES).
/// Invalid entries are silently skipped — check for ``None`` if needed.
///
///     mols = chematic.bulk.parse(["CCO", "c1ccccc1", "INVALID"])
///     valid = [m for m in mols if m is not None]
#[pyfunction]
pub fn parse(smiles: Vec<String>) -> Vec<Option<Mol>> {
    parse_parallel(&smiles)
        .into_iter()
        .map(|opt| {
            opt.map(|mol| Mol {
                inner: Arc::new(mol),
                props: Default::default(),
            })
        })
        .collect()
}

// ---------------------------------------------------------------------------
// bulk.ecfp4 — batch ECFP4 fingerprints as numpy (N, 2048) uint8
// ---------------------------------------------------------------------------

/// Compute ECFP4 fingerprints for a list of SMILES in parallel.
///
/// Returns a numpy array of shape ``(N, 2048)`` with ``dtype=uint8`` (0/1).
/// Invalid SMILES are silently skipped; ``N`` equals the number of valid molecules.
///
///     fps = chematic.bulk.ecfp4(smiles_list)  # shape (N, 2048)
///
///     # Use directly with scikit-learn:
///     from sklearn.decomposition import PCA
///     pca = PCA(n_components=50)
///     coords = pca.fit_transform(fps.astype(float))
#[pyfunction]
pub fn ecfp4<'py>(py: Python<'py>, smiles: Vec<String>) -> Bound<'py, PyArray2<u8>> {
    let bits: Vec<Vec<u8>> = smiles
        .par_iter()
        .filter_map(|s| chematic_smiles::parse(s).ok())
        .map(|mol| bitvec2048_to_bits(&chematic_fp::ecfp4(&mol)))
        .collect();

    let n = bits.len();
    if n == 0 {
        return Array2::<u8>::zeros((0, 2048)).into_pyarray(py);
    }
    let flat: Vec<u8> = bits.into_iter().flatten().collect();
    Array2::from_shape_vec((n, 2048), flat)
        .expect("shape mismatch")
        .into_pyarray(py)
}

// ---------------------------------------------------------------------------
// bulk.maccs — batch MACCS 166-bit as numpy (N, 166) uint8
// ---------------------------------------------------------------------------

/// Compute MACCS 166-bit fingerprints for a list of SMILES in parallel.
///
/// Returns a numpy array of shape ``(N, 166)`` with ``dtype=uint8`` (0/1).
#[pyfunction]
pub fn maccs<'py>(py: Python<'py>, smiles: Vec<String>) -> Bound<'py, PyArray2<u8>> {
    let bits: Vec<Vec<u8>> = smiles
        .par_iter()
        .filter_map(|s| chematic_smiles::parse(s).ok())
        .map(|mol| {
            let fp = chematic_fp::maccs(&mol);
            (0..166usize).map(|i| u8::from(fp.get(i))).collect()
        })
        .collect();

    let n = bits.len();
    if n == 0 {
        return Array2::<u8>::zeros((0, 166)).into_pyarray(py);
    }
    let flat: Vec<u8> = bits.into_iter().flatten().collect();
    Array2::from_shape_vec((n, 166), flat)
        .expect("shape mismatch")
        .into_pyarray(py)
}

// ---------------------------------------------------------------------------
// bulk.descriptors — parallel descriptor computation
// ---------------------------------------------------------------------------

/// Compute 55+ descriptors for a list of SMILES in parallel.
///
/// Returns a list of dicts. Each dict matches the output of ``mol.descriptors()``.
/// Invalid SMILES are silently skipped.
///
///     descs = chematic.bulk.descriptors(smiles_list)
///
///     # Build a Pandas DataFrame in one line:
///     import pandas as pd
///     df = pd.DataFrame(chematic.bulk.descriptors(smiles_list))
#[pyfunction]
pub fn descriptors<'py>(py: Python<'py>, smiles: Vec<String>) -> PyResult<Vec<Bound<'py, PyDict>>> {
    let records = smiles
        .par_iter()
        .filter_map(|smiles| chematic_smiles::parse(smiles).ok())
        .map(|molecule| compute_descriptor_record(&molecule))
        .collect::<Vec<_>>();
    records
        .into_iter()
        .map(|record| descriptor_record_to_pydict(py, record))
        .collect()
}

// ---------------------------------------------------------------------------
// bulk.descriptors_array — columnar numpy output
// ---------------------------------------------------------------------------

/// Compute descriptors and return selected columns as numpy arrays.
///
/// Returns a dict mapping each requested column name to a 1-D numpy array.
/// Float columns use ``float64``; bool columns use ``bool``; optional float
/// columns (``pka_acid``, ``pka_base``) use ``float64`` with ``NaN`` for None.
///
/// Raises ``ValueError`` for unknown column names.
///
///     result = chematic.bulk.descriptors_array(smiles, ["mw", "logp", "tpsa"])
///     df = pd.DataFrame(result)      # fast, no per-molecule dict allocation
///     mw = result["mw"]              # numpy.ndarray, dtype float64
#[pyfunction]
pub fn descriptors_array<'py>(
    py: Python<'py>,
    smiles: Vec<String>,
    columns: Vec<String>,
) -> PyResult<Bound<'py, pyo3::types::PyDict>> {
    let plan = DescriptorColumnPlan::new(&columns)?;
    let rows = smiles
        .par_iter()
        .filter_map(|smiles| chematic_smiles::parse(smiles).ok())
        .map(|molecule| SelectedDescriptorContext::new(&molecule, &plan).into_row())
        .collect::<Vec<_>>();
    materialize_descriptor_columns(py, &rows, &columns)
}

// ---------------------------------------------------------------------------
// bulk.tanimoto — pairwise Tanimoto similarity matrix
// ---------------------------------------------------------------------------

/// Compute an M×N Tanimoto similarity matrix (ECFP4) in parallel.
///
/// Args:
///     smiles_a: list of M SMILES strings (query molecules)
///     smiles_b: list of N SMILES strings (database molecules)
///
/// Returns:
///     numpy array of shape ``(M, N)`` with ``dtype=float32``.
///     Invalid SMILES are silently skipped.
///
///     # Self-similarity matrix:
///     sim = chematic.bulk.tanimoto(smiles, smiles)
///     # → (N, N) matrix; diagonal ≈ 1.0
///
///     # Query vs database:
///     sim = chematic.bulk.tanimoto(query_smiles, db_smiles)
///     # → (len(query_smiles), len(db_smiles))
#[pyfunction]
pub fn tanimoto<'py>(
    py: Python<'py>,
    smiles_a: Vec<String>,
    smiles_b: Vec<String>,
) -> Bound<'py, PyArray2<f32>> {
    let fps_a: Vec<chematic_fp::bitvec::BitVec2048> = smiles_a
        .par_iter()
        .filter_map(|s| chematic_smiles::parse(s).ok())
        .map(|mol| chematic_fp::ecfp4(&mol))
        .collect();

    let fps_b: Vec<chematic_fp::bitvec::BitVec2048> = smiles_b
        .par_iter()
        .filter_map(|s| chematic_smiles::parse(s).ok())
        .map(|mol| chematic_fp::ecfp4(&mol))
        .collect();

    let m = fps_a.len();
    let n = fps_b.len();

    if m == 0 || n == 0 {
        return Array2::<f32>::zeros((m, n)).into_pyarray(py);
    }

    // Re-use the existing Rust bulk Tanimoto matrix (row-major flat Vec)
    let flat = chematic_fp::bulk::tanimoto_matrix_parallel(&fps_a, &fps_b);
    Array2::from_shape_vec((m, n), flat)
        .expect("shape mismatch")
        .into_pyarray(py)
}

// ---------------------------------------------------------------------------
// bulk.tanimoto_search — 1:N similarity search
// ---------------------------------------------------------------------------

/// Compute Tanimoto similarity between one query molecule and N database molecules.
///
/// Returns a numpy array of shape ``(N,)`` with ``dtype=float32``.
///
///     scores = chematic.bulk.tanimoto_search("c1ccccc1", db_smiles)
///     top10_idx = np.argsort(scores)[::-1][:10]
#[pyfunction]
pub fn tanimoto_search<'py>(
    py: Python<'py>,
    query: String,
    smiles: Vec<String>,
) -> PyResult<Bound<'py, PyArray1<f32>>> {
    let query_mol = chematic_smiles::parse(&query)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;
    let query_fp = chematic_fp::ecfp4(&query_mol);

    let db_fps: Vec<chematic_fp::bitvec::BitVec2048> = smiles
        .par_iter()
        .filter_map(|s| chematic_smiles::parse(s).ok())
        .map(|mol| chematic_fp::ecfp4(&mol))
        .collect();

    if db_fps.is_empty() {
        return Ok(Array1::<f32>::zeros(0).into_pyarray(py));
    }

    let scores = chematic_fp::bulk::tanimoto_slice(&query_fp, &db_fps);
    Ok(Array1::from_vec(scores).into_pyarray(py))
}

// ---------------------------------------------------------------------------
// bulk.map4 — batch MAP4 fingerprints as numpy (N, 1024) uint32
// ---------------------------------------------------------------------------

/// Compute MAP4 fingerprints for a list of SMILES in parallel.
///
/// Returns a numpy array of shape ``(N, 1024)`` with ``dtype=uint32``.
/// Use :func:`chematic.tanimoto_map4` (not :func:`chematic.tanimoto`) for similarity.
/// Invalid SMILES are silently skipped.
///
///     fps = chematic.bulk.map4(smiles_list)  # shape (N, 1024), dtype uint32
#[pyfunction]
pub fn map4<'py>(py: Python<'py>, smiles: Vec<String>) -> Bound<'py, PyArray2<u32>> {
    let fps: Vec<Vec<u32>> = smiles
        .par_iter()
        .filter_map(|s| chematic_smiles::parse(s).ok())
        .map(|mol| chematic_fp::map4_default(&mol))
        .collect();

    let n = fps.len();
    if n == 0 {
        return Array2::<u32>::zeros((0, 1024)).into_pyarray(py);
    }
    let flat: Vec<u32> = fps.into_iter().flatten().collect();
    Array2::from_shape_vec((n, 1024), flat)
        .expect("shape mismatch")
        .into_pyarray(py)
}

// ---------------------------------------------------------------------------
// bulk.hdf — parallel Hyper-Dimensional Fingerprints
// ---------------------------------------------------------------------------

/// Compute HDF fingerprints for a list of SMILES in parallel.
///
/// Returns a numpy array of shape ``(N, dim)`` with ``dtype=float32`` where
/// each row is a unit-norm HD vector.  Invalid SMILES produce all-zero rows.
///
/// Args:
///     smiles: list of SMILES strings
///     dim: vector dimension (default 1024)
///     radius: neighborhood radius (default 2)
///     seed: reproducibility seed (default 42)
///
///     fps = chematic.bulk.hdf(["CCO", "c1ccccc1"], dim=512)
///     sims = fps @ fps.T          # cosine similarity matrix (N×N)
#[pyfunction]
#[pyo3(signature = (smiles, dim = 1024, radius = 2, seed = 42))]
pub fn hdf<'py>(
    py: Python<'py>,
    smiles: Vec<String>,
    dim: usize,
    radius: usize,
    seed: u64,
) -> Bound<'py, PyArray2<f32>> {
    let config = chematic_fp::HdfConfig { dim, radius, seed };
    let fps: Vec<Vec<f32>> = smiles
        .par_iter()
        .map(|s| {
            chematic_smiles::parse(s)
                .map(|mol| chematic_fp::hdf(&mol, &config).0)
                .unwrap_or_else(|_| vec![0.0f32; dim])
        })
        .collect();

    let n = fps.len();
    if n == 0 {
        return Array2::<f32>::zeros((0, dim)).into_pyarray(py);
    }
    let flat: Vec<f32> = fps.into_iter().flatten().collect();
    Array2::from_shape_vec((n, dim), flat)
        .expect("shape mismatch")
        .into_pyarray(py)
}

// ---------------------------------------------------------------------------
// bulk.substructure_search — parallel SMARTS screen
// ---------------------------------------------------------------------------

/// Screen a list of SMILES against a SMARTS pattern in parallel.
///
/// Returns a ``list[bool]`` of length N (same order as ``smiles``).
/// ``True`` if the molecule matches; ``False`` if it does not match or if the
/// SMILES is invalid.  The SMARTS is compiled once and shared across threads.
///
///     hits = chematic.bulk.substructure_search("[nH]", smiles_list)
///     actives = [s for s, h in zip(smiles_list, hits) if h]
#[pyfunction]
#[pyo3(signature = (smarts, smiles))]
pub fn substructure_search(smarts: &str, smiles: Vec<String>) -> PyResult<Vec<bool>> {
    let query = chematic_smarts::parse_smarts(smarts)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;
    // Stop at the first embedding per molecule instead of enumerating every
    // match — this is an existence screen, not a match-collection call.
    let config = chematic_smarts::MatchConfig {
        max_matches: Some(1),
        uniquify: false,
        ..chematic_smarts::MatchConfig::default()
    };

    let results: Vec<bool> = smiles
        .par_iter()
        .map(|smi| {
            chematic_smiles::parse(smi)
                .map(|mol| chematic_smarts::has_match_perceived(&query, &mol, &config))
                .unwrap_or(false)
        })
        .collect();

    Ok(results)
}

// ---------------------------------------------------------------------------
// bulk.substructure_match — parallel SMARTS screen over pre-parsed Mol objects
// ---------------------------------------------------------------------------

/// Screen a list of pre-parsed :class:`Mol` objects against a SMARTS pattern in parallel.
///
/// Unlike :func:`bulk.substructure_search` (which re-parses SMILES on every call),
/// this function accepts already-parsed molecules — avoiding re-parsing overhead
/// when screening the same collection multiple times.
///
/// Returns the **indices** of matching molecules in the input list.
///
///     mols = [m for m in chematic.bulk.parse(smiles_list) if m is not None]
///     hits = chematic.bulk.substructure_match("[OH]", mols)
///     # hits = [2, 5, 11, ...]  — indices into mols
#[pyfunction]
#[pyo3(signature = (smarts, mols))]
pub fn substructure_match(smarts: &str, mols: Vec<Mol>) -> PyResult<Vec<usize>> {
    let query = chematic_smarts::parse_smarts(smarts)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;
    // Stop at the first embedding per molecule instead of enumerating every
    // match — this is an existence screen, not a match-collection call.
    let config = chematic_smarts::MatchConfig {
        max_matches: Some(1),
        uniquify: false,
        ..chematic_smarts::MatchConfig::default()
    };

    let indices: Vec<usize> = mols
        .par_iter()
        .enumerate()
        .filter_map(|(i, m)| {
            if chematic_smarts::has_match_perceived(&query, &m.inner, &config) {
                Some(i)
            } else {
                None
            }
        })
        .collect();

    Ok(indices)
}

// ---------------------------------------------------------------------------
// bulk.generate_3d — parallel 3D coordinate generation
// ---------------------------------------------------------------------------

/// Generate 3D coordinates for a list of SMILES strings in parallel.
///
/// Returns a list of coordinate arrays ``[[x, y, z], ...]``, one per heavy
/// atom.  Invalid SMILES or molecules that fail 3D generation return ``None``.
///
/// Parameters
/// ----------
/// smiles : list[str]
///     Input SMILES strings.
/// method : {"etkdg", "dreiding"}
///     Force-field / algorithm used.  ``"etkdg"`` (default) applies the
///     ETKDG knowledge base (chair/envelope rings, 80 torsion rules).
///     ``"dreiding"`` is faster and uses the DREIDING force field only.
///
/// Examples
/// --------
///     coords = chematic.bulk.generate_3d(["CCO", "c1ccccc1"])
///     # coords[0] = [[x0,y0,z0], [x1,y1,z1], ...]  (one row per heavy atom)
#[pyfunction]
#[pyo3(signature = (smiles, *, method = "etkdg"))]
pub fn generate_3d(smiles: Vec<String>, method: &str) -> Vec<Option<Vec<[f64; 3]>>> {
    let use_etkdg = method != "dreiding";
    smiles
        .par_iter()
        .map(|s| {
            let mol = chematic_smiles::parse(s).ok()?;
            let coords = if use_etkdg {
                chematic_3d::generate_coords_etkdg(&mol)
            } else {
                chematic_3d::generate_and_minimize_dreiding(&mol)
            };
            let pts: Vec<[f64; 3]> = (0..mol.atom_count() as u32)
                .map(|i| {
                    let p = coords.get(chematic_core::AtomIdx(i));
                    [p.x, p.y, p.z]
                })
                .collect();
            Some(pts)
        })
        .collect()
}

// ---------------------------------------------------------------------------
// bulk.tanimoto_matrix — all-pairs Tanimoto similarity matrix
// ---------------------------------------------------------------------------

/// Compute the all-pairs ECFP4 Tanimoto similarity matrix for a list of SMILES.
///
/// Returns a numpy array of shape ``(N, N)`` with ``dtype=float32``.
/// Molecules that fail to parse are silently excluded; the matrix rows/columns
/// correspond only to successfully parsed molecules.
///
/// Parameters
/// ----------
/// smiles : list[str]
///     Input SMILES strings.
///
/// Examples
/// --------
///     mat = chematic.bulk.tanimoto_matrix(["CCO", "c1ccccc1", "CC(=O)O"])
///     # mat.shape == (3, 3);  mat[i, j] == Tanimoto(smiles[i], smiles[j])
#[pyfunction]
pub fn tanimoto_matrix<'py>(py: Python<'py>, smiles: Vec<String>) -> Bound<'py, PyArray2<f32>> {
    let fps: Vec<chematic_fp::bitvec::BitVec2048> = smiles
        .par_iter()
        .filter_map(|s| chematic_smiles::parse(s).ok())
        .map(|mol| chematic_fp::ecfp4(&mol))
        .collect();

    let n = fps.len();
    if n == 0 {
        return Array2::<f32>::zeros((0, 0)).into_pyarray(py);
    }

    let flat = chematic_fp::bulk::tanimoto_matrix_parallel(&fps, &fps);
    Array2::from_shape_vec((n, n), flat)
        .expect("shape mismatch in tanimoto_matrix")
        .into_pyarray(py)
}

// ---------------------------------------------------------------------------
// bulk.standardize — batch molecule standardization
// ---------------------------------------------------------------------------

/// Standardize a list of :class:`Mol` objects in parallel.
///
/// Applies the default standardization pipeline to each molecule:
/// largest-fragment selection, charge neutralization, and canonical tautomer.
/// All molecules succeed (standardization never fails); the output list has
/// the same length as the input.
///
/// Parameters
/// ----------
/// mols : list[Mol]
///     Pre-parsed molecules (e.g. from :func:`bulk.parse`).
///
/// Examples
/// --------
///     mols = chematic.bulk.parse(["CC(=O)[O-].[Na+]", "Oc1ccccc1"])
///     valid = [m for m in mols if m is not None]
///     std = chematic.bulk.standardize(valid)
///     # Na+ salt → neutral acid; phenol → canonical tautomer
#[pyfunction]
pub fn standardize(mols: Vec<Mol>) -> Vec<Mol> {
    let opts = chematic_chem::StandardizeOptions::default();
    mols.par_iter()
        .map(|m| {
            let s = chematic_chem::standardize(&m.inner, &opts);
            Mol::bare(s)
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Register the bulk submodule
// ---------------------------------------------------------------------------

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    register_pyfunctions!(
        m,
        parse,
        ecfp4,
        maccs,
        map4,
        hdf,
        descriptors,
        descriptors_array,
        tanimoto,
        tanimoto_search,
        substructure_search,
        substructure_match,
        generate_3d,
        tanimoto_matrix,
        standardize,
    )
}
