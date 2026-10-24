//! Command dispatch grouped by public CLI workflow.

use super::*;

pub(super) fn run(cli: Cli) -> Result<(), String> {
    match cli.command {
        Command::Fract {
            input,
            output,
            inspect,
            from_json,
            max_input_bytes,
        } => {
            let text = read_limited_input(input.as_ref(), max_input_bytes)?;
            let result = fract::convert(&text, inspect, from_json, max_input_bytes)?;
            write_output(output.as_ref(), &result)
        }
        Command::DlpolyHistory {
            input,
            output,
            inspect,
            from_json,
            max_input_bytes,
        } => {
            let text = read_limited_input(input.as_ref(), max_input_bytes)?;
            let record = if from_json {
                chematic_mol::dlpoly_history_from_json(&text)
            } else {
                chematic_mol::parse_dlpoly_history_with_limits(
                    &text,
                    &chematic_mol::DlpolyHistoryLimits {
                        max_input_bytes,
                        ..Default::default()
                    },
                )
            }
            .map_err(|e| e.to_string())?;
            let result = if inspect {
                serde_json::to_string_pretty(&record).map_err(|e| e.to_string())? + "\n"
            } else {
                chematic_mol::write_dlpoly_history(&record).map_err(|e| e.to_string())?
            };
            write_output(output.as_ref(), &result)
        }
        Command::NwchemGeometry {
            input,
            output,
            inspect,
            from_json,
            max_input_bytes,
        } => {
            let text = read_limited_input(input.as_ref(), max_input_bytes)?;
            let record = if from_json {
                chematic_mol::nwchem_geometry_from_json(&text)
            } else {
                chematic_mol::parse_nwchem_geometry_with_limits(
                    &text,
                    &chematic_mol::NwchemLimits {
                        max_input_bytes,
                        ..Default::default()
                    },
                )
            }
            .map_err(|e| e.to_string())?;
            let result = if inspect {
                serde_json::to_string_pretty(&record).map_err(|e| e.to_string())? + "\n"
            } else {
                chematic_mol::write_nwchem_geometry(&record).map_err(|e| e.to_string())?
            };
            write_output(output.as_ref(), &result)
        }
        Command::TinkerXyz {
            input,
            output,
            inspect,
            from_json,
            max_input_bytes,
        } => {
            let text = read_limited_input(input.as_ref(), max_input_bytes)?;
            let record = if from_json {
                chematic_mol::tinker_xyz_from_json(&text)
            } else {
                chematic_mol::parse_tinker_xyz_with_limits(
                    &text,
                    &chematic_mol::TinkerParseLimits {
                        max_input_bytes,
                        ..Default::default()
                    },
                )
            }
            .map_err(|e| e.to_string())?;
            let result = if inspect {
                serde_json::to_string(&record)
                    .map(|s| format!("{s}\n"))
                    .map_err(|e| e.to_string())?
            } else {
                chematic_mol::write_tinker_xyz(&record).map_err(|e| e.to_string())?
            };
            write_output(output.as_ref(), &result)
        }
        Command::Bgf {
            input,
            output,
            inspect,
            from_json,
            max_input_bytes,
        } => {
            let text = read_limited_input(input.as_ref(), max_input_bytes)?;
            let record = if from_json {
                chematic_mol::bgf_record_from_json(&text)
            } else {
                chematic_mol::parse_bgf_with_limits(
                    &text,
                    &chematic_mol::BgfLimits {
                        max_input_bytes,
                        ..Default::default()
                    },
                )
            }
            .map_err(|e| e.to_string())?;
            let result = if inspect {
                serde_json::to_string_pretty(&record).map_err(|e| e.to_string())? + "\n"
            } else {
                chematic_mol::write_bgf(&record).map_err(|e| e.to_string())?
            };
            write_output(output.as_ref(), &result)
        }
        Command::Molden {
            input,
            output,
            inspect,
            from_json,
            max_input_bytes,
        } => {
            let text = read_limited_input(input.as_ref(), max_input_bytes)?;
            let record = if from_json {
                chematic_mol::molden_record_from_json(&text)
            } else {
                chematic_mol::parse_molden_with_limits(
                    &text,
                    &chematic_mol::MoldenLimits {
                        max_input_bytes,
                        ..Default::default()
                    },
                )
            }
            .map_err(|e| e.to_string())?;
            let result = if inspect {
                serde_json::to_string(&record)
                    .map(|s| format!("{s}\n"))
                    .map_err(|e| e.to_string())?
            } else {
                chematic_mol::write_molden(&record).map_err(|e| e.to_string())?
            };
            write_output(output.as_ref(), &result)
        }
        Command::DlpolyConfig {
            input,
            output,
            inspect,
            from_json,
            max_input_bytes,
        } => {
            let text = read_limited_input(input.as_ref(), max_input_bytes)?;
            let record = if from_json {
                chematic_mol::dlpoly_config_from_json(&text)
            } else {
                chematic_mol::parse_dlpoly_config_with_limits(
                    &text,
                    &chematic_mol::DlpolyConfigLimits {
                        max_input_bytes,
                        ..Default::default()
                    },
                )
            }
            .map_err(|e| e.to_string())?;
            let result = if inspect {
                serde_json::to_string(&record)
                    .map(|s| format!("{s}\n"))
                    .map_err(|e| e.to_string())?
            } else {
                chematic_mol::write_dlpoly_config(&record).map_err(|e| e.to_string())?
            };
            write_output(output.as_ref(), &result)
        }
        Command::GjfDocument {
            input,
            output,
            inspect,
            max_input_bytes,
            max_jobs,
        } => {
            let text = read_limited_input(input.as_ref(), max_input_bytes)?;
            let result = conversion_format::gjf_document(&text, inspect, max_jobs)?;
            write_output(output.as_ref(), &result)
        }
        command @ (Command::ConvertWithReport { .. } | Command::Convert { .. }) => {
            run_conversion(command)
        }
        command @ (Command::Descriptors { .. }
        | Command::Parse { .. }
        | Command::Cxsmiles { .. }
        | Command::Fingerprint { .. }
        | Command::Similarity { .. }
        | Command::Substructure { .. }
        | Command::Standardize { .. }
        | Command::Report { .. }) => run_molecule(command),
        command @ (Command::Reaction { .. }
        | Command::ReactionMatch { .. }
        | Command::ReactionBalance { .. }
        | Command::ReactionFingerprint { .. }
        | Command::ReactionSimilarity { .. }) => run_reaction(command),
        command @ (Command::BatchReport { .. }
        | Command::BatchDescriptors { .. }
        | Command::BatchCanonicalize { .. }
        | Command::BatchFingerprints { .. }
        | Command::BatchStandardize { .. }
        | Command::BatchSimilarity { .. }
        | Command::BatchSubstructure { .. }
        | Command::BatchReactions { .. }) => run_batch(command),
    }
}

fn run_conversion(command: Command) -> Result<(), String> {
    match command {
        Command::ConvertWithReport {
            input_format,
            output_format,
            policy,
            input,
            output,
            max_input_bytes,
        } => {
            let text = read_limited_input(input.as_ref(), max_input_bytes)?;
            let mut limits = chematic_mol::ConversionLimits::default();
            limits.mol2.max_input_bytes = max_input_bytes;
            let json = chematic_mol::convert_with_report_json(
                &text,
                &input_format,
                &output_format,
                &policy,
                &limits,
            )
            .map_err(|error| error.to_string())?;
            let refused = serde_json::from_str::<serde_json::Value>(&json)
                .map_err(|error| error.to_string())?["status"]
                == "refused";
            write_output(output.as_ref(), &format!("{json}\n"))?;
            if refused {
                Err("conversion refused; see JSON report".into())
            } else {
                Ok(())
            }
        }
        Command::Convert {
            input_format,
            output_format,
            input,
            output,
            max_input_bytes,
        } => {
            let text = read_limited_input(input.as_ref(), max_input_bytes)?;
            let converted = convert_text(&text, &input_format, &output_format)?;
            write_output(output.as_ref(), &converted)
        }
        _ => unreachable!("conversion dispatcher received another command"),
    }
}

fn emit_json(json: String) -> Result<(), String> {
    write_output(None, &format!("{json}\n"))
}

fn run_molecule(command: Command) -> Result<(), String> {
    match command {
        Command::Descriptors { smiles } => emit_json(descriptors_json(&smiles)?),
        Command::Parse { smiles } => emit_json(parse_json(&smiles)?),
        Command::Cxsmiles { cxsmiles } => emit_json(cxsmiles_json(&cxsmiles)?),
        Command::Fingerprint { smiles, algorithm } => {
            emit_json(fingerprint_json(&smiles, &algorithm)?)
        }
        Command::Similarity {
            smiles_a,
            smiles_b,
            algorithm,
        } => emit_json(similarity_json(&smiles_a, &smiles_b, &algorithm)?),
        Command::Substructure {
            smiles,
            smarts,
            use_chirality,
            rdkit_parity,
        } => emit_json(substructure_json(
            &smiles,
            &smarts,
            use_chirality,
            rdkit_parity,
        )?),
        Command::Standardize {
            smiles,
            largest_fragment_only,
        } => emit_json(standardize_json(&smiles, largest_fragment_only)?),
        Command::Report { smiles } => emit_json(report_json(&smiles)?),
        _ => unreachable!("molecule dispatcher received another command"),
    }
}

fn run_reaction(command: Command) -> Result<(), String> {
    match command {
        Command::Reaction { reaction_smiles } => emit_json(reaction_json(&reaction_smiles)?),
        Command::ReactionMatch {
            reaction_smiles,
            query,
        } => emit_json(reaction_match_json(&reaction_smiles, &query)?),
        Command::ReactionBalance { reaction_smiles } => {
            emit_json(reaction_balance_json(&reaction_smiles)?)
        }
        Command::ReactionFingerprint {
            reaction_smiles,
            mode,
        } => emit_json(reaction_fingerprint_json(&reaction_smiles, &mode)?),
        Command::ReactionSimilarity {
            reaction_a,
            reaction_b,
        } => emit_json(reaction_similarity_json(&reaction_a, &reaction_b)?),
        _ => unreachable!("reaction dispatcher received another command"),
    }
}

fn read_batch(input: Option<&PathBuf>, limits: &BatchLimits) -> Result<String, String> {
    read_limited_input(input, limits.max_input_bytes)
}

fn run_batch(command: Command) -> Result<(), String> {
    match command {
        Command::BatchReport { input, limits } => {
            let text = read_batch(input.as_ref(), &limits)?;
            emit_json(batch_report_json(&text, &limits)?)
        }
        Command::BatchDescriptors { input, limits } => {
            let text = read_batch(input.as_ref(), &limits)?;
            emit_json(batch_descriptors_json(&text, &limits)?)
        }
        Command::BatchCanonicalize { input, limits } => {
            let text = read_batch(input.as_ref(), &limits)?;
            emit_json(batch_canonicalize_json(&text, &limits)?)
        }
        Command::BatchFingerprints {
            input,
            algorithm,
            limits,
        } => {
            let text = read_batch(input.as_ref(), &limits)?;
            emit_json(batch_fingerprints_json(&text, &algorithm, &limits)?)
        }
        Command::BatchStandardize {
            input,
            largest_fragment_only,
            limits,
        } => {
            let text = read_batch(input.as_ref(), &limits)?;
            emit_json(batch_standardize_json(
                &text,
                largest_fragment_only,
                &limits,
            )?)
        }
        Command::BatchSimilarity {
            input,
            algorithm,
            limits,
        } => {
            let text = read_batch(input.as_ref(), &limits)?;
            emit_json(batch_similarity_json(&text, &algorithm, &limits)?)
        }
        Command::BatchSubstructure {
            input,
            use_chirality,
            rdkit_parity,
            limits,
        } => {
            let text = read_batch(input.as_ref(), &limits)?;
            emit_json(batch_substructure_json(
                &text,
                &limits,
                use_chirality,
                rdkit_parity,
            )?)
        }
        Command::BatchReactions { input, limits } => {
            let text = read_batch(input.as_ref(), &limits)?;
            emit_json(batch_reactions_json(&text, &limits)?)
        }
        _ => unreachable!("batch dispatcher received another command"),
    }
}
