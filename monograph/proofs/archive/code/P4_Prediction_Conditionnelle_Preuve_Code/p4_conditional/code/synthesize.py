#!/usr/bin/env python3
"""Synthesize certificates, WITHOUT treating experimental success as proof."""
import json, struct
from pathlib import Path
from fractions import Fraction as F
from collections import deque
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'; OUT.mkdir(exist_ok=True)

def phi(x,b):
    return F(1,3)+x/(1+2*x) if b else F(1,3)+x/(3*(3-2*x))

def graph(p):
    D=1<<p; js=[D//2]; ids={js[0]:0}; edges=[]
    for j in js:
        row=[]
        for b in (0,1):
            # Exact nearest-grid rounding, half upward.
            z=D*phi(F(j,D),b)
            k=(2*z.numerator+z.denominator)//(2*z.denominator)
            if k not in ids: ids[k]=len(js); js.append(k)
            row.append(ids[k])
        edges.append(row)
    return js,edges

def main():
    points={F(1,2):''}; layer=[(F(1,2),'')]
    for depth in range(14):
        new=[]
        for x,w in layer:
            for b in (0,1):
                y=phi(x,b); v=w+str(b); new.append((y,v)); points.setdefault(y,v)
        layer=new
    sorted_pts=sorted(points.items())
    configs=[]
    for p in range(4,23):
        D=1<<p; eps=F(49,264*D); js,edges=graph(p)
        chosen=[]; last=None
        for x,w in sorted_pts:
            if last is None or x-last>8*eps:
                chosen.append((x,w)); last=x
        min_gap=min(b[0]-a[0] for a,b in zip(chosen,chosen[1:]))/4
        delta=min(eps/32,(min_gap-2*eps)/8,F(1,1000))
        cert={'version':1,'p':p,'epsilon':str(eps),'radius':str(4*eps),
              'states':js,'transitions':edges,'initial':0,
              'witness_words':[w for x,w in chosen],
              'lower_states':len(chosen),'upper_states':len(js),
              'lower_bits':(len(chosen)-1).bit_length(),
              'upper_bits':(len(js)-1).bit_length(),
              'robust_parameter_radius':str(delta),
              'minimum_witness_tv':str(min_gap),
              'witness_search_depth':14}
        (OUT/f'certificate_p{p:02d}.json').write_text(json.dumps(cert,indent=2)+'\n')
        # Runtime binary: 16 byte header + 8 bytes per immutable node + state bytes.
        B=cert['upper_bits']; nb=(B+7)//8
        payload=struct.pack('<4sBBHII',b'PCV1',p,nb,len(js),B,0)
        payload+=b''.join(struct.pack('<IHH',j,*e) for j,e in zip(js,edges))
        payload+=bytes(nb)
        (OUT/f'machine_p{p:02d}.bin').write_bytes(payload)
        configs.append({k:v for k,v in cert.items() if k not in ('states','transitions','witness_words')})
        print(p,len(js),len(chosen),cert['lower_bits'],B,flush=True)
    (OUT/'synthesis_summary.json').write_text(json.dumps({'configs':configs,'prefixes_enumerated':len(sorted_pts)},indent=2)+'\n')
if __name__=='__main__': main()
