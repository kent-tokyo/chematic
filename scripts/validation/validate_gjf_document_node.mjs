import {readFileSync} from 'node:fs';
import {pathToFileURL} from 'node:url';
const wasm = await import(pathToFileURL(process.argv[2]).href);
const records = [];
for (const c of JSON.parse(readFileSync(0,'utf8'))) {
  let result;
  try {
    const output = wasm.roundtrip_gjf_document(c.source);
    // Keep the raw JSON: JSON.stringify(JSON.parse(...)) changes -0 into 0.
    result = {ok:true,output,inspection:wasm.parse_gjf_document_json(c.source)};
  } catch (error) { result = {ok:false,error:String(error)}; }
  records.push({id:c.id,...result});
}
console.log(JSON.stringify(records));
