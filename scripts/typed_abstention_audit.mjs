#!/usr/bin/env node
/**
 * The probes of typed_abstention_audit.py on the WASM package: how each
 * kind of abstention (malformed / unsupported / ambiguous / resource_limit)
 * surfaces in JavaScript.
 *
 * Usage: node typed_abstention_audit.mjs --package-dir DIR --output OUT
 *   (DIR: an unpacked npm package; web target, initialised with initSync)
 */
import fs from "node:fs";
import path from "node:path";
import { fileURLToPath, pathToFileURL } from "node:url";

const arg = (name) => process.argv[process.argv.indexOf(name) + 1];
const dir = arg("--package-dir");
const wasm = await import(pathToFileURL(path.join(dir, "chematic_wasm.js")).href);
wasm.initSync({ module: fs.readFileSync(path.join(dir, "chematic_wasm_bg.wasm")) });
const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "..");
const nucleic = JSON.parse(fs.readFileSync(path.join(root, "validation/nucleic_acid_document_contract.json"), "utf8"));
const nucleicCase = (code) => nucleic.invalid_cases.find((c) => c.expected_code === code).document;
const overLimit = () => {
  const doc = { schema: "chematic.nucleic-acid.v1", atom_ids: [], strands: [], linkages: [], annotations: {} };
  for (let i = 0; i < 65; i += 1) {
    doc.atom_ids.push(`a${i}`);
    doc.strands.push({ id: `s${i}`, kind: "dna", annotations: {}, residues: [
      { id: `r${i}`, base: "A", sugar: "deoxyribose", atom_refs: [`a${i}`], annotations: {} }] });
  }
  return doc;
};
const pipeline = (smiles) => {
  const config = JSON.parse(wasm.pipeline_v2_stereo_safe_config_json("mmff94_bond_angle_strict", "diagnostic_only"));
  return wasm.embed_pipeline_v2_json(wasm.parse_smiles(smiles), JSON.stringify(config.config ?? config));
};

const probes = [
  ["smiles", "malformed", () => wasm.parse_smiles("C1CC(")],
  ["smarts", "malformed", () => wasm.match_smarts_smiles("CC", "[C;")],
  ["mol_block", "malformed", () => wasm.mol_block_coords_json("not a mol block")],
  ["smirks", "malformed", () => wasm.run_reactants_checked("[C:1>>", "CC", true)],
  ["smirks", "unsupported", () => wasm.run_reactants_checked("[C:1].[O:2]>>[C:1][O:2]", "C.O.N", true)],
  ["3d", "unsupported", () => pipeline("Cl[Pt](Cl)([NH3])[NH3]")],
  ["nucleic_acid", "ambiguous", () => wasm.nucleic_acid_validate_json(JSON.stringify(nucleicCase("ambiguous_atom_mapping")))],
  ["nucleic_acid", "unsupported", () => wasm.nucleic_acid_validate_json(JSON.stringify(nucleicCase("unsupported_linkage_topology")))],
  ["nucleic_acid", "resource_limit", () => wasm.nucleic_acid_validate_json(JSON.stringify(overLimit()))],
  ["mmcif", "resource_limit", () => wasm.mmcif_to_json("data_x\n".repeat(400000))],
];

function describe(value) {
  let parsed = value;
  if (typeof value === "string") {
    try { parsed = JSON.parse(value); } catch { return { outcome: "returned", fields: {} }; }
  }
  if (parsed && typeof parsed === "object" && !Array.isArray(parsed)) {
    const fields = {};
    for (const k of ["ok", "status", "reason", "code"]) if (k in parsed) fields[k] = parsed[k];
    if (parsed.error && typeof parsed.error === "object") fields["error.code"] = parsed.error.code ?? parsed.error.kind ?? parsed.error.cause?.kind;
    return { outcome: Object.keys(fields).length ? "envelope" : "returned", fields };
  }
  return { outcome: "returned", fields: {} };
}

const rows = probes.map(([family, category, call]) => {
  const row = { family, category };
  try {
    Object.assign(row, describe(call()));
  } catch (error) {
    const fields = {};
    for (const k of ["code", "category", "kind", "reason", "path"]) if (error && k in Object(error)) fields[k] = String(error[k]);
    Object.assign(row, { outcome: "raised", class: error?.constructor?.name ?? typeof error, fields,
      message: String(error?.message ?? error).slice(0, 160) });
  }
  return row;
});
const version = JSON.parse(fs.readFileSync(path.join(dir, "package.json"), "utf8")).version;
fs.writeFileSync(arg("--output"), JSON.stringify({ schema: "typed-abstention-audit/v1", binding: "wasm",
  chematic_version: version, probes: rows }, null, 1) + "\n");
for (const r of rows) console.log(r.family.padEnd(13), r.category.padEnd(15), r.outcome.padEnd(9), r.class ?? "", JSON.stringify(r.fields));
