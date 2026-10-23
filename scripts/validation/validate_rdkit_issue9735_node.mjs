import {readFileSync} from 'node:fs';
import {pathToFileURL} from 'node:url';
const wasm=await import(pathToFileURL(process.argv[2]).href);
const cases=JSON.parse(readFileSync(0,'utf8'));
const rows=[];
for(const c of cases){
  const mol=wasm.parse_smiles(c.smiles);let canonical;
  try{canonical={ok:true,value:mol.rdkit_smiles()};}
  catch(error){canonical={ok:false,error:String(error)};}
  finally{mol.free();}
  rows.push({id:c.id,canonical});
}
process.stdout.write(JSON.stringify({rows}));
