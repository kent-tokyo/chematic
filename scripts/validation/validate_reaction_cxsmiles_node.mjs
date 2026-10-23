import { createRequire } from 'node:module';
import { readFileSync, writeFileSync } from 'node:fs';
import { createHash } from 'node:crypto';
import { gzipSync } from 'node:zlib';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
const require=createRequire(import.meta.url);
const modulePath=path.resolve(process.argv[2]);
const outputPath=process.argv[3];
const wasm=require(modulePath);
const fixture=new URL('../../validation/rdkit-9732-2026-10-24/cases.json',import.meta.url);
const sha=p=>createHash('sha256').update(readFileSync(p)).digest('hex');
const available=typeof wasm.normalize_reaction_cxsmiles==='function';
const cases=[];
if(available) for(const c of JSON.parse(readFileSync(fixture)).positive) for(const canonical of [false,true]) {
  const row={id:c.id,input:c.input,canonical};
  try {row.output=wasm.normalize_reaction_cxsmiles(c.input,canonical);}
  catch(error) {row.error=String(error);}
  cases.push(row);
}
const result={schema:'chematic.reaction-cx-capture.v1',engine:'wasm',available,cases,
  runtime:{path:modulePath,sha256:sha(modulePath),wasm_sha256:sha(modulePath.replace(/\.js$/,'_bg.wasm'))},
  fixture_sha256:sha(fixture),observer_sha256:sha(fileURLToPath(import.meta.url)),utc:new Date().toISOString()};
writeFileSync(outputPath,gzipSync(JSON.stringify(result,null,2)+'\n'));
console.log(`WASM available: ${available}; captured: ${cases.length}`);
