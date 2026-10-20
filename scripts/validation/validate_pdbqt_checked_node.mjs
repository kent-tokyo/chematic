import fs from 'node:fs';
import {createRequire} from 'node:module';
const require=createRequire(import.meta.url);
const wasm=require(process.argv[2]);
const cases=JSON.parse(fs.readFileSync(0,'utf8'));
const rows=cases.map(c=>{
  try {
    return {id:c.id,output:wasm.smiles_to_pdbqt_checked(c.smiles,JSON.stringify(c.coords),JSON.stringify(c.charges),c.name)};
  } catch(e) {return {id:c.id,error:String(e)};}
});
process.stdout.write(JSON.stringify(rows));
