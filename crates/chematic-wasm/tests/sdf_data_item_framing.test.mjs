import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {pathToFileURL} from 'node:url';
const wasm = await import(process.env.CHEMATIC_WASM_MODULE
  ? pathToFileURL(process.env.CHEMATIC_WASM_MODULE).href
  : new URL('../pkg-node/chematic_wasm.js',import.meta.url).href);
const packet = JSON.parse(readFileSync(new URL('../../../validation/sdf_data_item_framing/cases.json',import.meta.url),'utf8'));
let checked = 0;
for (const ctab of Object.values(packet.ctabs)) {
  for (const ending of ['\n','\r\n']) {
    for (const entry of packet.cases) {
      const input = (ctab+entry.data+'$$$$\n').replaceAll('\n',ending);
      const records = JSON.parse(wasm.sdf_to_records_json(input));
      assert.equal(records.length,1);
      assert.deepEqual(records[0].properties,entry.expected_properties,entry.id);
      const batch = JSON.parse(wasm.sdf_records_batch_json(input,0,1));
      assert.equal(batch.record_count,1);
      assert.deepEqual(batch.records[0].record.properties,entry.expected_properties,entry.id);
      checked++;
    }
  }
}
console.log(`SDF data item framing: ${checked} V2000/V3000 LF/CRLF cases passed`);
