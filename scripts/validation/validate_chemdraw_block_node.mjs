import {readFileSync} from 'node:fs';
import {pathToFileURL} from 'node:url';
const wasm=await import(pathToFileURL(process.argv[2]).href);
const result=JSON.parse(readFileSync(0,'utf8')).map(c=>{
  const input=Buffer.from(c.input_hex,'hex');
  const report=JSON.parse(wasm.chemdraw_block_report_json(input,c.format,c.allow_metadata_loss));
  const blocks=wasm.mols_from_chemdraw_block(input,c.format,c.allow_metadata_loss).map(m=>{
    const block=wasm.to_mol_v3000_block(m);m.free();return block;
  });
  return {report,blocks};
});
process.stdout.write(JSON.stringify(result));
