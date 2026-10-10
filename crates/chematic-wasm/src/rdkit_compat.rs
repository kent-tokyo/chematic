//! JS-friendly bindings for the RDKit-compatible interoperability surface.
//!
//! Atom indices in every JSON result are zero-based and refer to the
//! [`MolHandle`] atom order. Coordinate arrays use that same order.

use std::collections::BTreeMap;
use std::rc::Rc;

use serde::Serialize;
use wasm_bindgen::prelude::*;

use crate::{MolHandle, bounded_json_string, enforce_wasm_input_len, enforce_wasm_molecule_size};

const MAX_3D_ATOMS: usize = 512;
const MAX_STEREOISOMERS: usize = 4096;
const MAX_ALIGNMENT_MATCHES: usize = 1_000_000;
const MAX_ALIGNMENT_ITERATIONS: u32 = 1_000_000;

fn js_error(code: &str, message: impl std::fmt::Display) -> JsValue {
    JsValue::from_str(&format!("{code}: {message}"))
}

fn chemistry_error(error: impl std::fmt::Display) -> JsValue {
    js_error("rdkit_compatibility_error", error)
}

fn validate_3d_size(mol: &MolHandle) -> Result<(), JsValue> {
    enforce_wasm_molecule_size(&mol.inner)?;
    if mol.inner.atom_count() > MAX_3D_ATOMS {
        return Err(js_error(
            "rdkit_3d_limit",
            format!(
                "molecule has {} atoms; the RDKit-compatible WASM 3D limit is {MAX_3D_ATOMS}",
                mol.inner.atom_count()
            ),
        ));
    }
    Ok(())
}

fn parse_json<T: serde::de::DeserializeOwned>(input: &str, label: &str) -> Result<T, JsValue> {
    enforce_wasm_input_len(label, input)?;
    serde_json::from_str(input).map_err(|error| {
        js_error(
            "rdkit_invalid_json",
            format!("{label} must be valid JSON: {error}"),
        )
    })
}

fn parse_coords(input: &str, label: &str, atom_count: usize) -> Result<Vec<[f64; 3]>, JsValue> {
    let coords: Vec<[f64; 3]> = parse_json(input, label)?;
    if coords.len() != atom_count {
        return Err(js_error(
            "rdkit_coordinate_count",
            format!(
                "{label} has {} rows, but the molecule has {atom_count} atoms",
                coords.len()
            ),
        ));
    }
    if coords.iter().flatten().any(|value| !value.is_finite()) {
        return Err(js_error(
            "rdkit_non_finite_coordinate",
            format!("{label} contains a non-finite coordinate"),
        ));
    }
    Ok(coords)
}

fn parse_optional_weights(
    input: Option<String>,
    expected_len: usize,
) -> Result<Option<Vec<f64>>, JsValue> {
    let Some(input) = input else {
        return Ok(None);
    };
    let weights: Vec<f64> = parse_json(&input, "weights_json")?;
    if weights.len() != expected_len {
        return Err(js_error(
            "rdkit_weight_count",
            format!(
                "weights_json has {} entries, expected {expected_len}",
                weights.len()
            ),
        ));
    }
    if weights
        .iter()
        .any(|weight| !weight.is_finite() || *weight <= 0.0)
    {
        return Err(js_error(
            "rdkit_invalid_weight",
            "weights must be finite and greater than zero",
        ));
    }
    Ok(Some(weights))
}

fn validate_max_matches(max_matches: usize) -> Result<(), JsValue> {
    if max_matches == 0 || max_matches > MAX_ALIGNMENT_MATCHES {
        return Err(js_error(
            "rdkit_alignment_limit",
            format!("max_matches must be between 1 and {MAX_ALIGNMENT_MATCHES}"),
        ));
    }
    Ok(())
}

fn validate_max_iterations(max_iterations: u32) -> Result<(), JsValue> {
    if max_iterations == 0 || max_iterations > MAX_ALIGNMENT_ITERATIONS {
        return Err(js_error(
            "rdkit_alignment_limit",
            format!("max_iterations must be between 1 and {MAX_ALIGNMENT_ITERATIONS}"),
        ));
    }
    Ok(())
}

#[derive(Serialize)]
struct AlignmentJson {
    rmsd: f64,
    transform: [[f64; 4]; 4],
    atom_map: Vec<(usize, usize)>,
}

/// Result of an RDKit-compatible PDB, XYZ, or MOL2 reader.
///
/// `coords_json()` is a JSON array with one `[x,y,z]` row per zero-based
/// molecule atom. `smiles` is empty for XYZ and unsanitized readers.
#[wasm_bindgen]
pub struct RdkitReadResult {
    molecule: Rc<chematic_core::Molecule>,
    coords: Vec<[f64; 3]>,
    smiles: String,
}

#[wasm_bindgen]
impl RdkitReadResult {
    /// A cloned handle to the parsed molecule.
    pub fn molecule(&self) -> MolHandle {
        MolHandle {
            inner: Rc::clone(&self.molecule),
        }
    }

    /// Coordinates as JSON, in zero-based molecule atom order.
    pub fn coords_json(&self) -> Result<String, JsValue> {
        bounded_json_string(&self.coords)
    }

    /// RDKit canonical SMILES when the reader sanitized the molecule.
    #[wasm_bindgen(getter)]
    pub fn smiles(&self) -> String {
        self.smiles.clone()
    }
}

#[wasm_bindgen]
impl MolHandle {
    /// RDKit-compatible SMARTS in input atom order.
    pub fn rdkit_smarts(
        &self,
        isomeric: bool,
        rooted_at_atom: Option<usize>,
    ) -> Result<String, JsValue> {
        enforce_wasm_molecule_size(&self.inner)?;
        chematic_smiles::rdkit_smarts(&self.inner, isomeric, rooted_at_atom)
            .map_err(chemistry_error)
    }

    /// RDKit-compatible CXSMARTS.
    pub fn rdkit_cx_smarts(&self) -> Result<String, JsValue> {
        enforce_wasm_molecule_size(&self.inner)?;
        chematic_smiles::rdkit_cx_smarts(&self.inner).map_err(chemistry_error)
    }

    /// RDKit-compatible PDB block. Optional coordinates are a JSON array in
    /// zero-based molecule atom order; omit them to write zero coordinates.
    pub fn rdkit_pdb_block(&self, coords_json: Option<String>) -> Result<String, JsValue> {
        enforce_wasm_molecule_size(&self.inner)?;
        let coords = coords_json
            .as_deref()
            .map(|json| parse_coords(json, "coords_json", self.inner.atom_count()))
            .transpose()?;
        chematic_smiles::rdkit_pdb_block(&self.inner, coords.as_deref()).map_err(chemistry_error)
    }

    /// RDKit-compatible Murcko scaffold SMILES.
    pub fn rdkit_murcko_scaffold(&self) -> Result<String, JsValue> {
        enforce_wasm_molecule_size(&self.inner)?;
        chematic_smiles::rdkit_murcko_scaffold(&self.inner).map_err(chemistry_error)
    }

    /// Chiral centres as `[[atomIndex,"R"|"S"|"?"], ...]` JSON.
    pub fn rdkit_chiral_centers_json(&self, include_unassigned: bool) -> Result<String, JsValue> {
        enforce_wasm_molecule_size(&self.inner)?;
        let centers = chematic_smiles::rdkit_chiral_centers(&self.inner, include_unassigned)
            .map_err(chemistry_error)?;
        bounded_json_string(&centers)
    }

    /// Number of default RDKit stereoisomers as a decimal string (u128-safe).
    pub fn rdkit_stereoisomer_count(&self) -> Result<String, JsValue> {
        enforce_wasm_molecule_size(&self.inner)?;
        chematic_smiles::rdkit_stereoisomer_count(&self.inner)
            .map(|count| count.to_string())
            .map_err(chemistry_error)
    }

    /// Sorted RDKit stereoisomer SMILES JSON. `max_isomers` is required and
    /// must be between 1 and 4096 so browser work remains bounded.
    pub fn rdkit_stereoisomer_smiles_json(&self, max_isomers: usize) -> Result<String, JsValue> {
        enforce_wasm_molecule_size(&self.inner)?;
        if max_isomers == 0 || max_isomers > MAX_STEREOISOMERS {
            return Err(js_error(
                "rdkit_stereoisomer_limit",
                format!("max_isomers must be between 1 and {MAX_STEREOISOMERS}"),
            ));
        }
        let smiles = chematic_smiles::rdkit_stereoisomer_smiles(&self.inner, max_isomers)
            .map_err(chemistry_error)?;
        bounded_json_string(&smiles)
    }

    /// RDKit MolHash by case-insensitive function name.
    pub fn rdkit_mol_hash(&self, function: &str, use_cx_smiles: bool) -> Result<String, JsValue> {
        enforce_wasm_molecule_size(&self.inner)?;
        let function =
            chematic_smiles::RdkitHashFunction::from_name(function).ok_or_else(|| {
                js_error(
                    "rdkit_hash_function",
                    format!("unknown RDKit hash function: {function}"),
                )
            })?;
        chematic_smiles::rdkit_mol_hash(&self.inner, function, use_cx_smiles)
            .map_err(chemistry_error)
    }

    /// RDKit Morgan folded bitInfo as `{bit:[[atomIndex,radius],...]}` JSON.
    pub fn rdkit_morgan_bit_info_json(
        &self,
        radius: u32,
        n_bits: usize,
        include_chirality: bool,
    ) -> Result<String, JsValue> {
        enforce_wasm_molecule_size(&self.inner)?;
        let config = rdkit_morgan_config(radius, n_bits, include_chirality)?;
        let result =
            chematic_fp::rdkit_morgan_fingerprint(&self.inner, &config).map_err(chemistry_error)?;
        let info: BTreeMap<usize, Vec<(u32, u32)>> = result
            .folded_bit_info
            .into_iter()
            .map(|(bit, mut environments)| {
                environments.sort_unstable();
                (bit, environments)
            })
            .collect();
        bounded_json_string(&info)
    }

    /// RDKit `AlignMol` result JSON. Coordinates and optional atom map and
    /// weights use zero-based molecule atom order.
    #[allow(clippy::too_many_arguments)]
    pub fn rdkit_align_json(
        &self,
        probe_coords_json: &str,
        reference_coords_json: &str,
        atom_map_json: Option<String>,
        weights_json: Option<String>,
        reflect: bool,
        max_iterations: u32,
    ) -> Result<String, JsValue> {
        validate_3d_size(self)?;
        validate_max_iterations(max_iterations)?;
        let probe = parse_coords(
            probe_coords_json,
            "probe_coords_json",
            self.inner.atom_count(),
        )?;
        let reference = parse_coords(
            reference_coords_json,
            "reference_coords_json",
            self.inner.atom_count(),
        )?;
        let atom_map: Option<Vec<(usize, usize)>> = atom_map_json
            .as_deref()
            .map(|json| parse_json(json, "atom_map_json"))
            .transpose()?;
        let weight_count = atom_map.as_ref().map_or(self.inner.atom_count(), Vec::len);
        let weights = parse_optional_weights(weights_json, weight_count)?;
        let result = chematic_smiles::rdkit_align_mol(
            &self.inner,
            &probe,
            &reference,
            atom_map.as_deref(),
            weights.as_deref(),
            reflect,
            max_iterations,
        )
        .map_err(chemistry_error)?;
        bounded_json_string(&AlignmentJson {
            rmsd: result.rmsd,
            transform: result.transform.rows(),
            atom_map: result.atom_map,
        })
    }

    /// RDKit `GetBestRMS`, bounded by `max_matches <= 1_000_000`.
    pub fn rdkit_best_rms(
        &self,
        probe_coords_json: &str,
        reference_coords_json: &str,
        max_matches: usize,
        symmetrize: bool,
        weights_json: Option<String>,
    ) -> Result<f64, JsValue> {
        validate_3d_size(self)?;
        validate_max_matches(max_matches)?;
        let probe = parse_coords(
            probe_coords_json,
            "probe_coords_json",
            self.inner.atom_count(),
        )?;
        let reference = parse_coords(
            reference_coords_json,
            "reference_coords_json",
            self.inner.atom_count(),
        )?;
        let weights = parse_optional_weights(weights_json, self.inner.atom_count())?;
        chematic_smiles::rdkit_best_rms(
            &self.inner,
            &probe,
            &reference,
            max_matches,
            symmetrize,
            weights.as_deref(),
        )
        .map(|result| result.rmsd)
        .map_err(chemistry_error)
    }

    /// RDKit `GetBestAlignmentTransform` as JSON, bounded by `max_matches`.
    pub fn rdkit_best_alignment_json(
        &self,
        probe_coords_json: &str,
        reference_coords_json: &str,
        max_matches: usize,
        symmetrize: bool,
        weights_json: Option<String>,
    ) -> Result<String, JsValue> {
        validate_3d_size(self)?;
        validate_max_matches(max_matches)?;
        let probe = parse_coords(
            probe_coords_json,
            "probe_coords_json",
            self.inner.atom_count(),
        )?;
        let reference = parse_coords(
            reference_coords_json,
            "reference_coords_json",
            self.inner.atom_count(),
        )?;
        let weights = parse_optional_weights(weights_json, self.inner.atom_count())?;
        let result = chematic_smiles::rdkit_best_rms(
            &self.inner,
            &probe,
            &reference,
            max_matches,
            symmetrize,
            weights.as_deref(),
        )
        .map_err(chemistry_error)?;
        bounded_json_string(&AlignmentJson {
            rmsd: result.rmsd,
            transform: result.transform.rows(),
            atom_map: result.atom_map,
        })
    }

    /// RDKit `CalcRMS`, bounded by `max_matches <= 1_000_000`.
    pub fn rdkit_calc_rms(
        &self,
        probe_coords_json: &str,
        reference_coords_json: &str,
        max_matches: usize,
        symmetrize: bool,
        weights_json: Option<String>,
    ) -> Result<f64, JsValue> {
        validate_3d_size(self)?;
        validate_max_matches(max_matches)?;
        let probe = parse_coords(
            probe_coords_json,
            "probe_coords_json",
            self.inner.atom_count(),
        )?;
        let reference = parse_coords(
            reference_coords_json,
            "reference_coords_json",
            self.inner.atom_count(),
        )?;
        let weights = parse_optional_weights(weights_json, self.inner.atom_count())?;
        chematic_smiles::rdkit_calc_rms(
            &self.inner,
            &probe,
            &reference,
            max_matches,
            symmetrize,
            weights.as_deref(),
        )
        .map_err(chemistry_error)
    }

    /// Seeded RDKit ETKDGv3 coordinates as JSON. Call on an
    /// explicit-hydrogen molecule. At most 512 atoms are accepted.
    pub fn rdkit_embed_json(
        &self,
        random_seed: i32,
        max_iterations: u32,
    ) -> Result<String, JsValue> {
        validate_3d_size(self)?;
        if max_iterations > MAX_ALIGNMENT_ITERATIONS {
            return Err(js_error(
                "rdkit_embed_limit",
                format!("max_iterations must not exceed {MAX_ALIGNMENT_ITERATIONS}"),
            ));
        }
        let options = chematic_3d::rdkit_embed::RdkitEmbedOptions {
            random_seed,
            max_iterations,
            ..Default::default()
        };
        let coords = chematic_3d::rdkit_embed::rdkit_embed_molecule(&self.inner, &options)
            .map_err(chemistry_error)?;
        bounded_json_string(&coords)
    }

    /// RDKit distance-geometry bounds matrix JSON. Row and column `i`
    /// refer to zero-based molecule atom `i`; upper bounds are above the
    /// diagonal and lower bounds below it. At most 512 atoms are accepted.
    pub fn rdkit_bounds_matrix_json(
        &self,
        set15bounds: bool,
        do_triangle_smoothing: bool,
        use_macrocycle_14_config: bool,
    ) -> Result<String, JsValue> {
        validate_3d_size(self)?;
        let matrix = chematic_3d::rdkit_embed::rdkit_bounds_matrix(
            &self.inner,
            set15bounds,
            do_triangle_smoothing,
            use_macrocycle_14_config,
        )
        .map_err(chemistry_error)?;
        bounded_json_string(&matrix)
    }
}

/// RDKit-compatible `MolFromPDBBlock` with molecule and coordinates.
#[wasm_bindgen]
pub fn rdkit_from_pdb_block(
    text: &str,
    sanitize: bool,
    remove_hs: bool,
    flavor: u32,
    proximity_bonding: bool,
) -> Result<RdkitReadResult, JsValue> {
    enforce_wasm_input_len("PDB input", text)?;
    let result = chematic_smiles::rdkit_mol_from_pdb_block(
        text,
        sanitize,
        remove_hs,
        flavor,
        proximity_bonding,
    )
    .map_err(chemistry_error)?
    .ok_or_else(|| {
        js_error(
            "rdkit_parse_failed",
            "RDKit returned no molecule for PDB input",
        )
    })?;
    enforce_wasm_molecule_size(&result.molecule)?;
    Ok(RdkitReadResult {
        molecule: Rc::new(result.molecule),
        coords: result.coords,
        smiles: result.smiles,
    })
}

/// RDKit-compatible `MolFromXYZBlock`; the molecule has atoms but no bonds.
#[wasm_bindgen]
pub fn rdkit_from_xyz_block(text: &str) -> Result<RdkitReadResult, JsValue> {
    enforce_wasm_input_len("XYZ input", text)?;
    let (molecule, coords) =
        chematic_smiles::rdkit_mol_from_xyz_block(text).map_err(chemistry_error)?;
    enforce_wasm_molecule_size(&molecule)?;
    Ok(RdkitReadResult {
        molecule: Rc::new(molecule),
        coords,
        smiles: String::new(),
    })
}

/// RDKit-compatible `MolFromMol2Block` with molecule and coordinates.
#[wasm_bindgen]
pub fn rdkit_from_mol2_block(
    text: &str,
    sanitize: bool,
    remove_hs: bool,
    cleanup_substructures: bool,
) -> Result<RdkitReadResult, JsValue> {
    enforce_wasm_input_len("MOL2 input", text)?;
    let result = chematic_smiles::rdkit_mol_from_mol2_block(
        text,
        sanitize,
        remove_hs,
        cleanup_substructures,
    )
    .map_err(chemistry_error)?;
    enforce_wasm_molecule_size(&result.molecule)?;
    Ok(RdkitReadResult {
        molecule: Rc::new(result.molecule),
        coords: result.coords,
        smiles: result.smiles,
    })
}

fn rdkit_morgan_config(
    radius: u32,
    n_bits: usize,
    include_chirality: bool,
) -> Result<chematic_fp::RdkitMorganConfig, JsValue> {
    let radius = match radius {
        0 => chematic_fp::RdkitMorganRadius::R0,
        1 => chematic_fp::RdkitMorganRadius::R1,
        2 => chematic_fp::RdkitMorganRadius::R2,
        3 => chematic_fp::RdkitMorganRadius::R3,
        _ => {
            return Err(js_error(
                "rdkit_morgan_config",
                "radius must be 0, 1, 2, or 3",
            ));
        }
    };
    let fp_size = match n_bits {
        128 => chematic_fp::RdkitMorganFpSize::B128,
        256 => chematic_fp::RdkitMorganFpSize::B256,
        512 => chematic_fp::RdkitMorganFpSize::B512,
        1024 => chematic_fp::RdkitMorganFpSize::B1024,
        2048 => chematic_fp::RdkitMorganFpSize::B2048,
        _ => {
            return Err(js_error(
                "rdkit_morgan_config",
                "n_bits must be 128, 256, 512, 1024, or 2048",
            ));
        }
    };
    Ok(chematic_fp::RdkitMorganConfig {
        radius,
        fp_size,
        include_chirality,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn molecule(smiles: &str) -> MolHandle {
        MolHandle {
            inner: Rc::new(chematic_smiles::parse(smiles).unwrap()),
        }
    }

    #[test]
    fn rdkit_writer_and_structure_json_match_rust_api() {
        let mol = molecule("O=C1CC[C@H](C)N1CCc1ccccc1");
        assert_eq!(
            mol.rdkit_smarts(true, None).unwrap(),
            chematic_smiles::rdkit_smarts(&mol.inner, true, None).unwrap()
        );
        assert_eq!(mol.rdkit_murcko_scaffold().unwrap(), "O=C1CCCN1CCc1ccccc1");
        let centers: Vec<(usize, String)> =
            serde_json::from_str(&mol.rdkit_chiral_centers_json(true).unwrap()).unwrap();
        assert_eq!(centers, vec![(4, "S".into())]);
    }

    #[test]
    fn rdkit_reader_results_preserve_atom_order_and_coords() {
        let result = rdkit_from_xyz_block("2\nethane\nC 0 0 0\nC 1.5 0 0\n").unwrap();
        assert_eq!(result.molecule().atom_count(), 2);
        let coords: Vec<[f64; 3]> = serde_json::from_str(&result.coords_json().unwrap()).unwrap();
        assert_eq!(coords, vec![[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]]);
    }

    #[test]
    fn rdkit_alignment_matches_rust_api() {
        let mol = molecule("CCO");
        let probe = "[[1,0,0],[2,0,0],[3,0,0]]";
        let reference = "[[0,0,0],[1,0,0],[2,0,0]]";
        let json = mol
            .rdkit_best_alignment_json(probe, reference, 100, true, None)
            .unwrap();
        let value: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(value["rmsd"].as_f64().unwrap() < 1e-12);
        assert_eq!(value["atom_map"].as_array().unwrap().len(), 3);
    }

    #[test]
    fn rdkit_embed_and_bounds_have_atom_aligned_shapes() {
        let mol = molecule("[H]OC([H])([H])C([H])([H])[H]");
        let coords: Vec<[f64; 3]> =
            serde_json::from_str(&mol.rdkit_embed_json(42, 0).unwrap()).unwrap();
        let bounds: Vec<Vec<f64>> =
            serde_json::from_str(&mol.rdkit_bounds_matrix_json(true, true, false).unwrap())
                .unwrap();
        assert_eq!(coords.len(), mol.atom_count());
        assert_eq!(bounds.len(), mol.atom_count());
        assert!(bounds.iter().all(|row| row.len() == mol.atom_count()));
    }
}
