import {readFileSync} from 'node:fs';
import {pathToFileURL} from 'node:url';
const wasm=await import(pathToFileURL(process.argv[2]).href);
const out=[];
for(const c of JSON.parse(readFileSync(0,'utf8'))){
  const row={id:c.id};
  for(const api of c.expected==='write'?['checked','legacy']:['checked']){
    try{row[api]={ok:true,output:(api==='checked'?wasm.write_pqr_json_checked:wasm.write_pqr_json)(c.records_json)};}
    catch(error){row[api]={ok:false,error:String(error)};}
  }
  out.push(row);
}
console.log(JSON.stringify(out));
