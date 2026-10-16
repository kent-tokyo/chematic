// Batch driver for validate_pdb_conect.py; handles and frames are both checked.
import { readFileSync } from 'node:fs';
import { pathToFileURL } from 'node:url';
const wasm = await import(pathToFileURL(process.argv[2]).href);
const cases = JSON.parse(readFileSync(0, 'utf8'));
const result = cases.map(c => {
  const frame = JSON.parse(wasm.pdb_conect_frame_json(c.input, c.extended));
  const handle = wasm.mol_from_pdb_conect(c.input, c.extended);
  const mol_block = wasm.to_mol_v3000_block(handle);
  handle.free();
  return {frame, mol_block};
});
process.stdout.write(JSON.stringify(result));
