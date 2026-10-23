import {readFileSync} from 'node:fs';
import {pathToFileURL} from 'node:url';
const wasm=await import(pathToFileURL(process.argv[2]).href);
const rows=[];
for(const c of JSON.parse(readFileSync(0,'utf8'))){
  const mol=wasm.parse_smiles(c.smiles),row={id:c.id};
  try{
    const before=wasm.to_moljson(mol);
    row.labelled={ok:true,formula:mol.formula_with_isotopes()};
    row.legacy={ok:true,formula:mol.formula()};
    row.graph_unchanged=wasm.to_moljson(mol)===before;
  }finally{mol.free();}
  rows.push(row);
}
process.stdout.write(JSON.stringify({rows}));
