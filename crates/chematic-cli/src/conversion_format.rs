/// Normalize the CLI's aliases for supported molecule formats.
pub(super) fn format_name(format: &str) -> Option<String> {
    match format
        .to_ascii_lowercase()
        .trim_start_matches('.')
        .replace('-', "_")
        .as_str()
    {
        "smi" | "smiles" => Some("smiles".to_string()),
        "mol" | "sdf" => Some("mol".to_string()),
        "v3000" | "mol_v3000" => Some("mol_v3000".to_string()),
        "gjf" | "com" => Some("gjf".to_string()),
        "mol2" => Some("mol2".to_string()),
        "cml" => Some("cml".to_string()),
        "cjson" => Some("cjson".to_string()),
        "moljson" => Some("moljson".to_string()),
        "cdxml" => Some("cdxml".to_string()),
        _ => None,
    }
}

/// Handle checked source projections before entering the generic graph bridge.
pub(super) fn source_format_conversion(
    text: &str,
    input: &str,
    output: &str,
) -> Option<Result<String, String>> {
    if input == "cjson" && output == "cml" {
        return Some(
            chematic_mol::project_cjson_to_cml_with_limits(
                text,
                &chematic_mol::ConversionLimits {
                    mol2: chematic_mol::Mol2ParseLimits {
                        max_atoms: super::MAX_SMILES_ATOMS,
                        ..Default::default()
                    },
                    ..Default::default()
                },
            )
            .map_err(|error| error.to_string()),
        );
    }
    if input != "gjf" && output != "gjf" {
        return None;
    }
    if input != output {
        return Some(Err(
            "GJF/COM conversion supports same-format source retention only".into(),
        ));
    }
    Some(
        chematic_mol::parse_gjf_source_record_with_limits(
            text,
            &chematic_mol::GaussianParseLimits {
                max_atoms: super::MAX_SMILES_ATOMS,
                ..Default::default()
            },
        )
        .map(|record| chematic_mol::write_gjf_source_record(&record))
        .map_err(|error| error.to_string()),
    )
}

pub(super) fn gjf_document(text: &str, inspect: bool, max_jobs: usize) -> Result<String, String> {
    if max_jobs > 128 {
        return Err("GJF document job limit must not exceed 128".into());
    }
    let document = chematic_mol::parse_gjf_source_document_with_limits(
        text,
        &chematic_mol::GjfDocumentLimits {
            gaussian: chematic_mol::GaussianParseLimits {
                max_atoms: super::MAX_SMILES_ATOMS,
                ..Default::default()
            },
            max_jobs,
        },
    )
    .map_err(|error| error.to_string())?;
    Ok(if inspect {
        format!("{}\n", document.to_json())
    } else {
        chematic_mol::write_gjf_source_document(&document)
    })
}
