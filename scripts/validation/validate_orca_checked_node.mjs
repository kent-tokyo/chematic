import {readFileSync} from 'node:fs';
import {pathToFileURL} from 'node:url';
const wasm=await import(pathToFileURL(process.argv[2]).href);
const cases=JSON.parse(readFileSync(0,'utf8'));
const records=[];
for(const c of cases){
  const record={id:c.id};
  const calls={checked:()=>wasm.write_orca_input_json_checked(c.input_json)};
  if(c.expected==='accept')calls.legacy=()=>wasm.write_orca_input_json(c.input_json.replace('"kind": "xyz"','"type": "xyz"'));
  for(const [name,call] of Object.entries(calls)){
    try{record[name]={ok:true,output:call()};}
    catch(error){record[name]={ok:false,error:String(error)};}
  }
  records.push(record);
}
process.stdout.write(JSON.stringify(records));
