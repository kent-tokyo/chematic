import {readFileSync} from 'node:fs';
import {pathToFileURL} from 'node:url';
const wasm=await import(pathToFileURL(process.argv[2]).href);
const records=[];
for(const c of JSON.parse(readFileSync(0,'utf8'))){
  const row={id:c.id};
  const functions={core:wasm.orca_input_to_json,legacy_coords:wasm.orca_input_coords_json};
  if(process.argv[3]==='candidate')functions.checked=wasm.orca_cartesian_geometry_json;
  for(const [name,call] of Object.entries(functions)){
    try{row[name]={ok:true,value:call(c.source)};}
    catch(error){row[name]={ok:false,error:String(error)};}
  }
  records.push(row);
}
// Keep returned JSON strings intact so -0 is not normalized by JSON.stringify.
console.log(JSON.stringify(records));
