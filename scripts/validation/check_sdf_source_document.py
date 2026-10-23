#!/usr/bin/env python3
"""Bounded development comparison of complete SDF retention and record accounting."""
from __future__ import annotations
import argparse,hashlib,importlib,json,tempfile,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from scripts.checks.check_interchange_preservation import bounded_run,identity,require,strict_json
NODE="""
const fs=require('node:fs'),w=require(process.argv[1]);
const cases=JSON.parse(fs.readFileSync(0,'utf8'));
process.stdout.write(JSON.stringify(cases.map(c=>{try{return {id:c.id,output:w.convert_common_format(c.input,'sdf','sdf')}}catch(e){return {id:c.id,error:String(e)}}})));
"""
def framing(text):
    """Limited independent physical framing observer, not a CTAB/chemistry parser."""
    records=[];pending=[];unexpected=False;complete=True
    def consume(lines):
        nonlocal unexpected,complete
        if not '\n'.join(lines).strip():return
        rows=lines
        if len(rows)<4 or not rows[3].rstrip().endswith(('V2000','V3000')):
            unexpected=True;return
        records.append(rows[0]);complete=complete and any(row=='M  END' for row in rows[4:])
    for raw in text.split('\n'):
        line=raw.removesuffix('\r')
        if line=='$$$$':consume(pending);pending=[]
        else:pending.append(line)
    consume(pending)
    return dict(record_count=len(records),titles=records,complete_terminators=complete,unexpected_text=unexpected)

def observe(source,output,titles):
    if output is None:return dict(exact_source=False,record_count=False,ordered_titles=False,clean_framing=False)
    parsed=framing(output)
    return dict(exact_source=source==output,record_count=parsed['record_count']==len(titles),ordered_titles=parsed['titles']==titles,clean_framing=bool(parsed['record_count']) and parsed['complete_terminators'] and not parsed['unexpected_text'])
def terminal(raw):return {k:v.decode('utf-8',errors='replace') if isinstance(v,bytes) else v for k,v in raw.items()}
def framed_carrier(text):
    newline='\r\n' if '\r\n' in text else '\n'
    if not text.endswith('\n'):text+=newline
    if text.split('\n')[-2].removesuffix('\r')!='$$$$':
        if not text.endswith(newline+newline):text+=newline
        text+='$$$$'+newline
    return text

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--wasm-module',type=Path,required=True);parser.add_argument('--obabel',type=Path,default=Path('/opt/homebrew/bin/obabel'));parser.add_argument('--output',type=Path,required=True);parser.add_argument('--baseline',action='store_true');args=parser.parse_args()
    import chematic
    extension=importlib.import_module('chematic.chematic');packet_path=ROOT/'validation/sdf_source_document/cases.json';packet=strict_json(packet_path.read_text(encoding='utf-8'));cases=[]
    for c in packet['cases']:
        require(hashlib.sha256(c['input'].encode()).hexdigest()==c['sha256'],'fixture hash mismatch')
        require(framing(c['input'])['titles']==c['titles'],'fixture framing mismatch')
        for layout,text in [('lf',c['input']),('crlf',c['input'].replace('\n','\r\n')),('no-final-newline',c['input'].rstrip('\n'))]:cases.append(dict(id=c['id']+'/'+layout,input=text,titles=c['titles']))
    raw=bounded_run(['node','-e',NODE,str(args.wasm_module.resolve())],json.dumps(cases).encode());require(raw['resource_status']=='completed' and raw['exit_code']==0,'Node failed');node=strict_json(raw['stdout']);require(len(node)==len(cases),'Node count mismatch');results=[]
    for c,w in zip(cases,node):
        require(c['id']==w['id'],'Node order mismatch');row=dict(id=c['id'],input_sha256=hashlib.sha256(c['input'].encode()).hexdigest(),expected_titles=c['titles'],wasm=w)
        try:row['python']=dict(output=chematic.convert_format(c['input'],'sdf','sdf'))
        except ValueError as error:row['python']=dict(error=str(error))
        for label,options in [('openbabel',[]),('openbabel_v3000',['-x3'])]:
            command=[str(args.obabel.resolve()),'-isdf','-osdf']+options;r=terminal(bounded_run(command,c['input'].encode()));r['command']=command
            if r['resource_status']=='completed' and r['exit_code']==0:r['output']=r['stdout']
            row[label]=r
        carrier=framed_carrier(c['input'])
        with tempfile.TemporaryDirectory(prefix='sdf-copy-') as directory:
            infile=Path(directory)/'source.sdf';infile.write_bytes(carrier.encode())
            for label,options in [('openbabel_copy',['-aT']),('openbabel_copy_full',[])]:
                command=[str(args.obabel.resolve()),'-isdf',str(infile),'-ocopy']+options;r=terminal(bounded_run(command));r['command']=command;r['input_sha256']=hashlib.sha256(carrier.encode()).hexdigest();r['workflow']='Seekable explicitly SDF-framed raw-copy carrier; title-only read' if options else 'Seekable explicitly SDF-framed raw-copy carrier; normal molecular read'
                if r['resource_status']=='completed' and r['exit_code']==0:r['output']=r['stdout']
                r['observation']=observe(carrier,r.get('output'),c['titles']);row[label]=r
        for label in ['python','wasm','openbabel','openbabel_v3000']:
            row[label]['observation']=observe(c['input'],row[label].get('output'),c['titles'])
            if row[label].get('output') is not None:row[label]['framing']=framing(row[label]['output'])
        results.append(row)
    labels=['python','wasm','openbabel','openbabel_v3000','openbabel_copy','openbabel_copy_full'];summary={label:{metric:sum(r[label]['observation'][metric] for r in results) for metric in results[0][label]['observation']} for label in labels}
    record=dict(schema='chematic.sdf-source-document-results.v1',baseline=args.baseline,packet=identity(packet_path),python_package=identity(chematic.__file__),python_extension=identity(extension.__file__),wasm_module=identity(args.wasm_module),wasm_binary=identity(args.wasm_module.with_name('chematic_wasm_bg.wasm')),openbabel=identity(args.obabel),openbabel_version=terminal(bounded_run([str(args.obabel.resolve()),'-V'])),case_count=len(results),summary=summary,cases=results)
    args.output.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8');print(json.dumps(summary))
    if not args.baseline:require(all(all(r[label]['observation'].values()) for r in results for label in labels[:2]),'whole SDF retention failed')
if __name__=='__main__':main()
