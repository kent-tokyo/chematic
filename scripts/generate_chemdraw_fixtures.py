#!/usr/bin/env python3
"""Original concrete CDX fixtures, independently paired with authored CDXML."""
import json
import math
import struct
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]

def prop(tag,data):return struct.pack('<HH',tag,len(data))+data

def obj(tag,ident,body):return struct.pack('<HI',tag,ident)+body+b'\0\0'

def encode(fragments,header=28,metadata=False,reverse=False):
    parts=[];xml=['<CDXML BondLength="30"><page id="2">']
    for fi,(atoms,bonds) in enumerate(fragments):
        body=[];xml.append(f'<fragment id="{100+fi}">')
        for i,(element,x,y,charge,isotope,hs) in enumerate(atoms):
            fields=[prop(0x0200,struct.pack('<ii',round(y*65536),round(x*65536))),prop(0x0402,struct.pack('<H',element))]
            attrs=f'id="{10+i}" Element="{element}" p="{x} {y}"'
            if charge:
                fields.append(prop(0x0421,struct.pack('<b',charge)));attrs+=f' Charge="{charge}"'
            if isotope:
                fields.append(prop(0x0420,struct.pack('<H',isotope)));attrs+=f' Isotope="{isotope}"'
            if hs is not None:
                fields.append(prop(0x042b,struct.pack('<H',hs)));attrs+=f' NumHydrogens="{hs}"'
            if reverse:fields.reverse()
            body.append(obj(0x8004,10+i,b''.join(fields)));xml.append(f'<n {attrs}/>')
        for i,(a,b,order) in enumerate(bonds):
            fields=[prop(0x0604,struct.pack('<I',10+a)),prop(0x0605,struct.pack('<I',10+b)),prop(0x0600,struct.pack('<H',{1:1,2:2,3:4,1.5:128}[order]))]
            if reverse:fields.reverse()
            body.append(obj(0x8005,1000+i,b''.join(fields)));xml.append(f'<b id="{1000+i}" B="{10+a}" E="{10+b}" Order="{order}"/>')
        parts.append(obj(0x8003,100+fi,b''.join(body)));xml.append('</fragment>')
    xml.append('</page></CDXML>')
    extra=prop(0x0805,struct.pack('<i',30*65536)) if metadata else b''
    data=b'VjCD0100\x04\x03\x02\x01'+b'\0'*(header-12)+obj(0x8000,1,extra+obj(0x8001,2,b''.join(parts)))+b'\0\0'
    return data,'\n'.join(xml)

def cases():
    def a(z,x,y=0,charge=0,isotope=0,hs=None):return(z,x,y,charge,isotope,hs)
    ring=[a(6,round(30*math.cos(i*math.pi/3)*16)/16,round(30*math.sin(i*math.pi/3)*16)/16,hs=1) for i in range(6)]
    examples=[('ethanol',[([a(6,0,hs=3),a(6,30,hs=2),a(8,60,hs=1)],[(0,1,1),(1,2,1)])]),
      ('ammonium',[([a(7,-1.25,2.5,1,0,4)],[])]),
      ('isotope-triple',[([a(6,0,0,0,13,3),a(6,30,hs=0),a(7,60,hs=0)],[(0,1,1),(1,2,3)])]),
      ('acetate',[([a(8,0,0,-1,0,0),a(6,30,hs=0),a(8,30,30,0,0,0),a(6,60,hs=3)],[(0,1,1),(1,2,2),(1,3,1)])]),
      ('aromatic',[(ring,[(i,(i+1)%6,1.5) for i in range(6)])]),
      ('alkene',[([a(9,0,30,hs=0),a(6,30,hs=1),a(6,60,hs=1),a(9,90,-30,hs=0)],[(0,1,1),(1,2,2),(2,3,1)])]),
      ('multi-fragment',[([a(11,0,0,1,0,0)],[]),([a(17,60,0,-1,0,0)],[])])]
    for name,fragments in examples:
        for header in (22,28):
            for reverse in (False,True):
                for metadata in (False,True):
                    cdx,xml=encode(fragments,header,metadata,reverse)
                    yield dict(id=f'{name}-{header}-{reverse}-{metadata}',cdx_hex=cdx.hex(),cdxml=xml,allow_metadata_loss=metadata)

if __name__=='__main__':
    folder=ROOT/'validation/rdkit_issues/chemdraw-block';folder.mkdir(parents=True,exist_ok=True)
    (folder/'cases.json').write_text(json.dumps(dict(schema='chematic.chemdraw-fixtures.v1',cases=list(cases())),indent=2)+'\n')
