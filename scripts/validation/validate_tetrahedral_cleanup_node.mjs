import fs from 'node:fs';
import {createRequire} from 'node:module';
const require=createRequire(import.meta.url);const w=require(process.argv[2]);
const cases=JSON.parse(fs.readFileSync(0,'utf8'));
const rows=cases.map(c=>{
 const m=w.parse_smiles(c.smiles);const r=w.remove_redundant_tetrahedral_stereo(m,c.atom_index,1000000);const out=r.molecule;
 const row={id:c.id,report:JSON.parse(r.report_json()),source_smiles:m.canonical_smiles(),output_smiles:out.canonical_smiles()};
 out.free();r.free();m.free();return row;
});
process.stdout.write(JSON.stringify(rows));
