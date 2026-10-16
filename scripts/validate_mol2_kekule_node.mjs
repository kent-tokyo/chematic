import { readFileSync } from 'node:fs';
import { createRequire } from 'node:module';
const wasm = createRequire(import.meta.url)(process.argv[2]);
const rows = JSON.parse(readFileSync(0, 'utf8')).map(text => {
  try { return {output: wasm.kekulize_mol2_record(text)}; }
  catch (error) { return {error: String(error)}; }
});
process.stdout.write(JSON.stringify(rows));
