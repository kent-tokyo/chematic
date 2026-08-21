// Plain node:assert unit tests for the Local Compound Explorer's pure
// functions -- same lightweight convention as
// crates/chematic-wasm/tests/*.test.mjs (no test framework).
//
// Run: node demo/explorer/tests/explorer.test.mjs

import assert from "node:assert/strict";
import path from "node:path";
import { fileURLToPath } from "node:url";

const __dirname = path.dirname(fileURLToPath(import.meta.url));

const { parseCsvText, detectColumns, csvRowsToRawRecords, parseSmiFileText } =
  await import(path.join(__dirname, "..", "parser.js"));
const { applyFilters, buildComparator, matchesFreeText } =
  await import(path.join(__dirname, "..", "table.js"));
const { csvField, exportToCsv, CSV_COLUMNS } =
  await import(path.join(__dirname, "..", "export.js"));

// ---------------------------------------------------------------------------
// parser.js: parseCsvText
// ---------------------------------------------------------------------------

console.log("=== parseCsvText ===");
{
  const rows = parseCsvText("name,smiles\naspirin,CC(=O)O\n");
  assert.deepEqual(rows, [["name", "smiles"], ["aspirin", "CC(=O)O"]]);
}
{
  // quoted field containing a comma
  const rows = parseCsvText('name,smiles\n"a, b",CCO\n');
  assert.deepEqual(rows, [["name", "smiles"], ["a, b", "CCO"]]);
}
{
  // quoted field containing an embedded real newline
  const rows = parseCsvText('name,notes\n"line1\nline2",x\n');
  assert.deepEqual(rows, [["name", "notes"], ["line1\nline2", "x"]]);
}
{
  // doubled-quote escaping
  const rows = parseCsvText('name,notes\n"He said ""hi""",x\n');
  assert.deepEqual(rows, [["name", "notes"], ['He said "hi"', "x"]]);
}
{
  // no trailing newline
  const rows = parseCsvText("a,b\n1,2");
  assert.deepEqual(rows, [["a", "b"], ["1", "2"]]);
}
{
  // CRLF line endings
  const rows = parseCsvText("a,b\r\n1,2\r\n");
  assert.deepEqual(rows, [["a", "b"], ["1", "2"]]);
}
{
  // bare CR line endings
  const rows = parseCsvText("a,b\r1,2\r");
  assert.deepEqual(rows, [["a", "b"], ["1", "2"]]);
}
{
  // BOM stripped
  const rows = parseCsvText("﻿a,b\n1,2\n");
  assert.deepEqual(rows, [["a", "b"], ["1", "2"]]);
}
console.log("parseCsvText: all assertions passed");

// ---------------------------------------------------------------------------
// parser.js: detectColumns
// ---------------------------------------------------------------------------

console.log("=== detectColumns ===");
{
  assert.deepEqual(detectColumns(["name", "smiles"]), { smilesCol: 1, nameCol: 0 });
  assert.deepEqual(detectColumns(["ID", "SMILES"]), { smilesCol: 1, nameCol: 0 });
  assert.deepEqual(detectColumns(["compound", "canonical_smiles"]), { smilesCol: 1, nameCol: 0 });
  assert.deepEqual(detectColumns(["Structure", "Name"]), { smilesCol: 0, nameCol: 1 });
  assert.deepEqual(detectColumns(["foo", "bar"]), { smilesCol: null, nameCol: null });
  // case-insensitive
  assert.deepEqual(detectColumns(["SmIlEs"]), { smilesCol: 0, nameCol: null });
}
console.log("detectColumns: all assertions passed");

// ---------------------------------------------------------------------------
// parser.js: csvRowsToRawRecords
// ---------------------------------------------------------------------------

console.log("=== csvRowsToRawRecords ===");
{
  const rows = [["name", "smiles"], ["aspirin", "CC(=O)O"], ["", ""], ["", "CCO"]];
  const records = csvRowsToRawRecords(rows, 1, 0);
  assert.deepEqual(records, [
    { name: "aspirin", smiles: "CC(=O)O" },
    { name: "", smiles: "CCO" },
  ]);
}
{
  // no name column
  const rows = [["smiles"], ["CCO"]];
  const records = csvRowsToRawRecords(rows, 0, null);
  assert.deepEqual(records, [{ name: "", smiles: "CCO" }]);
}
console.log("csvRowsToRawRecords: all assertions passed");

// ---------------------------------------------------------------------------
// parser.js: parseSmiFileText
// ---------------------------------------------------------------------------

console.log("=== parseSmiFileText ===");
{
  const records = parseSmiFileText("CCO\tethanol\nc1ccccc1  benzene\n\nCC(=O)O\n");
  assert.deepEqual(records, [
    { smiles: "CCO", name: "ethanol" },
    { smiles: "c1ccccc1", name: "benzene" },
    { smiles: "CC(=O)O", name: "" },
  ]);
}
console.log("parseSmiFileText: all assertions passed");

// ---------------------------------------------------------------------------
// table.js: matchesFreeText / applyFilters / buildComparator
// ---------------------------------------------------------------------------

console.log("=== matchesFreeText ===");
{
  const rec = { name: "Aspirin", inputSmiles: "CC(=O)O", canonicalSmiles: "CC(=O)O", formula: "C9H8O4" };
  assert.equal(matchesFreeText(rec, "aspirin"), true);
  assert.equal(matchesFreeText(rec, "ASPIRIN"), true);
  assert.equal(matchesFreeText(rec, "c9h8o4"), true);
  assert.equal(matchesFreeText(rec, "ibuprofen"), false);
  assert.equal(matchesFreeText(rec, ""), true);
  assert.equal(matchesFreeText(rec, "   "), true);
}
console.log("matchesFreeText: all assertions passed");

function makeRecord(overrides) {
  return {
    index: 0, name: "mol", inputSmiles: "C", canonicalSmiles: "C", formula: "CH4",
    status: "ok",
    descriptors: {
      mw: 100, logP: 1, tpsa: 20, hbd: 1, hba: 1, rotatableBonds: 0, qed: 0.5,
      lipinskiPasses: true, painsPasses: true, veberPasses: true, eganPasses: true,
      ghosePasses: true, reosPasses: true,
    },
    painsAlerts: [], similarity: null, errorMessage: null,
    ...overrides,
  };
}

console.log("=== applyFilters ===");
{
  const records = [
    makeRecord({ index: 0, name: "A", descriptors: { ...makeRecord().descriptors, mw: 100, logP: 1, lipinskiPasses: true, painsPasses: true } }),
    makeRecord({ index: 1, name: "B", descriptors: { ...makeRecord().descriptors, mw: 500, logP: 6, lipinskiPasses: false, painsPasses: false } }),
    makeRecord({ index: 2, name: "C", status: "error", descriptors: null, errorMessage: "bad smiles" }),
  ];

  // no filters -> everything passes
  assert.equal(applyFilters(records, {}).length, 3);

  // validOnly excludes the error row
  assert.deepEqual(applyFilters(records, { validOnly: true }).map((r) => r.index), [0, 1]);

  // lipinskiPass
  assert.deepEqual(applyFilters(records, { lipinskiPass: true }).map((r) => r.index), [0]);

  // painsPass
  assert.deepEqual(applyFilters(records, { painsPass: true }).map((r) => r.index), [0]);

  // mw range
  assert.deepEqual(applyFilters(records, { mwMin: 200 }).map((r) => r.index), [1]);
  assert.deepEqual(applyFilters(records, { mwMax: 200 }).map((r) => r.index), [0]);

  // combined filters
  assert.deepEqual(
    applyFilters(records, { validOnly: true, mwMax: 200, lipinskiPass: true }).map((r) => r.index),
    [0]
  );

  // free text
  assert.deepEqual(applyFilters(records, { text: "B" }).map((r) => r.index), [1]);

  // similarityMin excludes records with null similarity
  const withSim = [
    makeRecord({ index: 0, similarity: 0.9 }),
    makeRecord({ index: 1, similarity: 0.1 }),
    makeRecord({ index: 2, similarity: null }),
  ];
  assert.deepEqual(applyFilters(withSim, { similarityMin: 0.5 }).map((r) => r.index), [0]);
}
console.log("applyFilters: all assertions passed");

console.log("=== buildComparator ===");
{
  const records = [
    makeRecord({ index: 0, name: "Charlie", descriptors: { ...makeRecord().descriptors, mw: 300 } }),
    makeRecord({ index: 1, name: "Alice", descriptors: { ...makeRecord().descriptors, mw: 100 } }),
    makeRecord({ index: 2, name: "Bob", descriptors: { ...makeRecord().descriptors, mw: 200 } }),
  ];

  const byNameAsc = [...records].sort(buildComparator("name", "asc"));
  assert.deepEqual(byNameAsc.map((r) => r.name), ["Alice", "Bob", "Charlie"]);

  const byNameDesc = [...records].sort(buildComparator("name", "desc"));
  assert.deepEqual(byNameDesc.map((r) => r.name), ["Charlie", "Bob", "Alice"]);

  const byMwAsc = [...records].sort(buildComparator("mw", "asc"));
  assert.deepEqual(byMwAsc.map((r) => r.mw ?? r.descriptors.mw), [100, 200, 300]);

  // similarity: nulls always sort last, regardless of direction
  const withSim = [
    makeRecord({ index: 0, name: "A", similarity: null }),
    makeRecord({ index: 1, name: "B", similarity: 0.2 }),
    makeRecord({ index: 2, name: "C", similarity: 0.8 }),
  ];
  assert.deepEqual(
    [...withSim].sort(buildComparator("similarity", "desc")).map((r) => r.name),
    ["C", "B", "A"]
  );
  assert.deepEqual(
    [...withSim].sort(buildComparator("similarity", "asc")).map((r) => r.name),
    ["B", "C", "A"]
  );

  // inputOrder
  const byInputOrder = [...records].sort(buildComparator("inputOrder", "asc"));
  assert.deepEqual(byInputOrder.map((r) => r.index), [0, 1, 2]);
}
console.log("buildComparator: all assertions passed");

// ---------------------------------------------------------------------------
// export.js: csvField / exportToCsv
// ---------------------------------------------------------------------------

console.log("=== csvField ===");
{
  // a legitimate negative number must NOT be guarded (not in GUARDED_COLUMNS)
  assert.equal(csvField(-1.03, "logp"), "-1.03");
  assert.equal(csvField(-8.6477, "mw"), "-8.6477");

  // a guarded string column starting with '=' gets the OWASP apostrophe guard
  assert.equal(csvField("=SUM(A1)", "name"), "'=SUM(A1)");
  assert.equal(csvField("+1234", "input_smiles"), "'+1234");
  assert.equal(csvField("-1234", "formula"), "'-1234");
  assert.equal(csvField("@cmd", "error"), "'@cmd");
  assert.equal(csvField("\tx", "parse_status"), "'\tx");

  // ordinary strings are untouched
  assert.equal(csvField("aspirin", "name"), "aspirin");

  // RFC4180 quoting applies independently of the guard
  assert.equal(csvField("a,b", "name"), '"a,b"');
  assert.equal(csvField('say "hi"', "name"), '"say ""hi"""');
  assert.equal(csvField("=SUM(A1),B1", "name"), '"\'=SUM(A1),B1"');

  // null/undefined become empty string
  assert.equal(csvField(null, "similarity"), "");
  assert.equal(csvField(undefined, "similarity"), "");

  // booleans render as their string form, no guard
  assert.equal(csvField(true, "lipinski_passes"), "true");
  assert.equal(csvField(false, "pains_passes"), "false");
}
console.log("csvField: all assertions passed");

console.log("=== exportToCsv ===");
{
  const records = [
    makeRecord({
      index: 0, name: "Caffeine", inputSmiles: "Cn1cnc2c1c(=O)n(C)c(=O)n2C", canonicalSmiles: "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
      formula: "C8H10N4O2",
      descriptors: { mw: 194.19, logP: -1.03, tpsa: 61.8, hbd: 0, hba: 3, rotatableBonds: 0, qed: 0.6,
        lipinskiPasses: true, painsPasses: true, veberPasses: true, eganPasses: true, ghosePasses: true, reosPasses: true },
      similarity: 0.42,
    }),
    makeRecord({ index: 1, name: "Bad", status: "error", descriptors: null, errorMessage: "unexpected char" }),
  ];

  const csv = exportToCsv(records);
  const lines = csv.trim().split("\r\n");
  assert.equal(lines.length, 3); // header + 2 records
  assert.equal(lines[0], CSV_COLUMNS.join(","));

  const row0 = lines[1].split(",");
  assert.equal(row0[0], "0"); // input_index
  assert.equal(row0[1], "Caffeine"); // name
  assert.equal(row0[5], "194.19"); // mw
  assert.equal(row0[6], "-1.03"); // logp -- unguarded negative number
  assert.equal(row0[14], "0.42"); // similarity

  const row1 = lines[2].split(",");
  assert.equal(row1[1], "Bad");
  assert.equal(row1[14], ""); // similarity empty when null
  assert.equal(row1[15], "error"); // parse_status
  assert.equal(row1[16], "unexpected char"); // error
}
console.log("exportToCsv: all assertions passed");

console.log("\nexplorer.test.mjs: all assertions passed");
