import {readFileSync} from 'node:fs';
import {pathToFileURL} from 'node:url';
const wasm=await import(pathToFileURL(process.argv[2]).href);
const cases=JSON.parse(readFileSync(0,'utf8'));
const result=[];
for(const c of cases){
  const record={id:c.id};
  for(const [name,call] of [['direct',()=>wasm.roundtrip_gjf_record(c.source)],['common',()=>wasm.convert_common_format(c.source,'.COM','gjf')]]){
    try{record[name]={ok:true,output:call()};}
    catch(error){record[name]={ok:false,error:String(error)};}
  }
  result.push(record);
}
console.log(JSON.stringify(result));
