use serde_json::{json, Value};
use chematic_smiles::parse;
use chematic_chem::*;
use chematic_perception::*;

pub struct ChematicTools;

impl ChematicTools {
    pub fn new() -> Self {
        ChematicTools
    }

    pub async fn call_tool(&self, name: &str, args: &Value) -> Result<Value, String> {
        match name {
            "parse_smiles" => self.parse_smiles(args).await,
            "molecule_descriptors" => self.molecule_descriptors(args).await,
            "drug_likeness" => self.drug_likeness(args).await,
            "substructure_search" => self.substructure_search(args).await,
            "similarity" => self.similarity(args).await,
            "functional_groups" => self.functional_groups(args).await,
            "run_reaction" => self.run_reaction(args).await,
            "generate_3d" => self.generate_3d(args).await,
            "mcs" => self.mcs(args).await,
            "canonical_smiles" => self.canonical_smiles(args).await,
            "fingerprint" => self.fingerprint(args).await,
            "screen_library" => self.screen_library(args).await,
            _ => Err(format!("Unknown tool: {}", name)),
        }
    }

    async fn parse_smiles(&self, args: &Value) -> Result<Value, String> {
        let smiles = args["smiles"]
            .as_str()
            .ok_or("Missing SMILES argument")?;

        let mol = parse(smiles).map_err(|e| e.to_string())?;
        let canonical = chematic_smiles::canonical_smiles(&mol);
        let mw = molecular_weight(&mol);
        let formula = chematic_smiles::write(&mol);

        Ok(json!({
            "valid": true,
            "canonical": canonical,
            "molecular_weight": mw,
            "atom_count": mol.atom_count(),
            "bond_count": mol.bond_count(),
        }))
    }

    async fn molecule_descriptors(&self, args: &Value) -> Result<Value, String> {
        let smiles = args["smiles"]
            .as_str()
            .ok_or("Missing SMILES argument")?;

        let mol = parse(smiles).map_err(|e| e.to_string())?;

        Ok(json!({
            "molecular_weight": molecular_weight(&mol),
            "tpsa": tpsa(&mol),
            "logp": logp_crippen(&mol),
            "qed": qed(&mol),
            "hbd": hbd_count(&mol),
            "hba": hba_count(&mol),
            "heavy_atoms": heavy_atom_count(&mol),
            "rotatable_bonds": rotatable_bond_count(&mol),
            "aromatic_rings": aromatic_ring_count(&mol),
        }))
    }

    async fn drug_likeness(&self, args: &Value) -> Result<Value, String> {
        let smiles = args["smiles"]
            .as_str()
            .ok_or("Missing SMILES argument")?;

        let mol = parse(smiles).map_err(|e| e.to_string())?;

        Ok(json!({
            "lipinski": lipinski_passes(&mol),
            "veber": veber_passes(&mol),
            "egan": egan_passes(&mol),
            "reos": reos_passes(&mol),
            "ghose": ghose_passes(&mol),
        }))
    }

    async fn substructure_search(&self, args: &Value) -> Result<Value, String> {
        let smarts = args["smarts"]
            .as_str()
            .ok_or("Missing SMARTS argument")?;
        let smiles = args["smiles"]
            .as_str()
            .ok_or("Missing SMILES argument")?;

        let mol = parse(smiles).map_err(|e| e.to_string())?;
        let query = chematic_smarts::parse(smarts).map_err(|e| e.to_string())?;
        let matches = chematic_smarts::find_matches(&mol, &query);

        let match_arrays: Vec<Vec<u32>> = matches
            .iter()
            .map(|m| m.iter().map(|a| a.0).collect())
            .collect();

        Ok(json!({
            "matches": match_arrays,
            "match_count": match_arrays.len(),
        }))
    }

    async fn similarity(&self, args: &Value) -> Result<Value, String> {
        let smiles1 = args["smiles1"]
            .as_str()
            .ok_or("Missing smiles1 argument")?;
        let smiles2 = args["smiles2"]
            .as_str()
            .ok_or("Missing smiles2 argument")?;

        let mol1 = parse(smiles1).map_err(|e| e.to_string())?;
        let mol2 = parse(smiles2).map_err(|e| e.to_string())?;

        let tanimoto = chematic_fp::tanimoto_ecfp4(&mol1, &mol2);

        Ok(json!({
            "tanimoto": tanimoto,
        }))
    }

    async fn functional_groups(&self, args: &Value) -> Result<Value, String> {
        let smiles = args["smiles"]
            .as_str()
            .ok_or("Missing SMILES argument")?;

        let mol = parse(smiles).map_err(|e| e.to_string())?;
        let groups = detect_functional_groups(&mol);

        let group_json: Vec<Value> = groups
            .iter()
            .map(|fg| {
                json!({
                    "name": fg.name,
                    "atom_count": fg.atoms.len(),
                })
            })
            .collect();

        Ok(json!({
            "functional_groups": group_json,
            "count": group_json.len(),
        }))
    }

    async fn run_reaction(&self, args: &Value) -> Result<Value, String> {
        let smirks = args["smirks"]
            .as_str()
            .ok_or("Missing SMIRKS argument")?;
        let smiles = args["smiles"]
            .as_str()
            .ok_or("Missing SMILES argument")?;

        let reactant = parse(smiles).map_err(|e| e.to_string())?;
        let rxn = chematic_rxn::parse_smirks(smirks).map_err(|e| e.to_string())?;
        let products = chematic_rxn::run_reactants(&reactant, &rxn);

        let product_smiles: Vec<Vec<String>> = products
            .iter()
            .map(|product_set| {
                product_set
                    .iter()
                    .map(|mol| chematic_smiles::write(mol))
                    .collect()
            })
            .collect();

        Ok(json!({
            "products": product_smiles,
        }))
    }

    async fn generate_3d(&self, args: &Value) -> Result<Value, String> {
        let smiles = args["smiles"]
            .as_str()
            .ok_or("Missing SMILES argument")?;

        let mol = parse(smiles).map_err(|e| e.to_string())?;
        let coords = chematic_3d::generate_3d_minimized(&mol, Default::default())
            .map_err(|e| e.to_string())?;

        let pdb = chematic_mol::to_pdb_block(&mol, &coords);

        Ok(json!({
            "pdb": pdb,
        }))
    }

    async fn mcs(&self, args: &Value) -> Result<Value, String> {
        let smiles_list = args["smiles_list"]
            .as_array()
            .ok_or("Missing smiles_list argument")?;

        let molecules: Result<Vec<_>, _> = smiles_list
            .iter()
            .map(|s| {
                s.as_str()
                    .ok_or("Invalid SMILES in list")
                    .and_then(|smi| parse(smi).map_err(|e| e.to_string()))
            })
            .collect();

        let mols = molecules?;
        if mols.is_empty() {
            return Err("Empty molecule list".to_string());
        }

        let mcs = chematic_smarts::find_mcs(&mols);

        Ok(json!({
            "mcs_smarts": mcs.to_string(),
            "atom_count": mcs.atom_count(),
        }))
    }

    async fn canonical_smiles(&self, args: &Value) -> Result<Value, String> {
        let smiles = args["smiles"]
            .as_str()
            .ok_or("Missing SMILES argument")?;

        let mol = parse(smiles).map_err(|e| e.to_string())?;
        let canonical = chematic_smiles::canonical_smiles(&mol);
        let inchi = chematic_inchi::inchi(&mol);
        let inchikey = chematic_inchi::inchikey(&mol);

        Ok(json!({
            "canonical": canonical,
            "inchi": inchi,
            "inchikey": inchikey,
        }))
    }

    async fn fingerprint(&self, args: &Value) -> Result<Value, String> {
        let smiles = args["smiles"]
            .as_str()
            .ok_or("Missing SMILES argument")?;
        let fp_type = args["fp_type"]
            .as_str()
            .unwrap_or("ecfp4");

        let mol = parse(smiles).map_err(|e| e.to_string())?;

        match fp_type {
            "ecfp4" => {
                let fp = chematic_fp::ecfp4_counts(&mol);
                let json_fp = serde_json::to_value(fp).map_err(|e| e.to_string())?;
                Ok(json!({
                    "type": "ecfp4",
                    "fingerprint": json_fp,
                }))
            }
            "maccs" => {
                let fp = chematic_fp::maccs_bitvec(&mol);
                let hex = format!("{:0width$x}", fp, width = 256 / 4);
                Ok(json!({
                    "type": "maccs",
                    "fingerprint_hex": hex,
                    "bit_count": fp.count_ones(),
                }))
            }
            "topo" => {
                let fp = chematic_fp::topo_path_bitvec(&mol);
                let hex = format!("{:x}", fp);
                Ok(json!({
                    "type": "topo",
                    "fingerprint_hex": hex,
                    "bit_count": fp.count_ones(),
                }))
            }
            _ => Err(format!("Unknown fingerprint type: {}", fp_type)),
        }
    }

    async fn screen_library(&self, args: &Value) -> Result<Value, String> {
        let smiles_list = args["smiles_list"]
            .as_array()
            .ok_or("Missing smiles_list argument")?;
        let smarts = args["smarts"]
            .as_str()
            .ok_or("Missing SMARTS argument")?;

        let query = chematic_smarts::parse(smarts).map_err(|e| e.to_string())?;

        let mut hits = Vec::new();
        for smi_val in smiles_list {
            let smi = smi_val
                .as_str()
                .ok_or("Invalid SMILES in list")?;
            let mol = parse(smi).map_err(|e| e.to_string())?;
            let matches = chematic_smarts::find_matches(&mol, &query);
            if !matches.is_empty() {
                hits.push(smi.to_string());
            }
        }

        Ok(json!({
            "hits": hits,
            "hit_count": hits.len(),
        }))
    }
}
