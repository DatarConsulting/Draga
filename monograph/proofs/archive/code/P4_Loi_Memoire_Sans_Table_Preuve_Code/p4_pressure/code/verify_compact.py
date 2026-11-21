#!/usr/bin/env python3
from pathlib import Path
from fractions import Fraction as F
import json,ctypes as C,subprocess
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results'
M=(((F(1,2),F(1,4)),(F(1,12),F(1,6))),((F(1,6),F(1,12)),(F(1,4),F(1,2))))
J=(F(3,8),F(5,8))
def bayes(x,b):
 a,b0=M[b][0];c,d=M[b][1]
 v=(1-x)*b0+x*d;u=(1-x)*a+x*c
 return v/(u+v)
def value(s,x):
 # s is newest first; chronological updates run in reverse order.
 for b in reversed(s):x=bayes(x,int(b))
 return x
def bounds(s):return value(s,J[0]),value(s,J[1])
def get_nodes(p):
 eps=F(49,132*(1<<p)); nodes=[]
 def rec(s):
  lo,hi=bounds(s);isleaf=hi-lo<=4*eps
  nodes.append({'word':s,'lo':lo,'hi':hi,'leaf':isleaf})
  if not isleaf:rec(s+'0');rec(s+'1')
 rec('');return nodes

def main():
 subprocess.run(['gcc','-std=c11','-O2','-Wall','-Wextra','-Werror','-Wno-misleading-indentation','-shared','-fPIC',str(ROOT/'code/compact_filter.c'),'-o',str(ROOT/'code/libcompact.so')],check=True)
 lib=C.CDLL(str(ROOT/'code/libcompact.so'));U8=C.POINTER(C.c_uint8);U64=C.POINTER(C.c_uint64)
 lib.cf_count.argtypes=[C.c_uint,U64,C.POINTER(C.c_uint)]
 for variant in ['rank','path']:
  getattr(lib,f'cf_{variant}_step').argtypes=[C.c_uint,U8,C.c_uint,C.c_uint]
  getattr(lib,f'cf_{variant}_output').argtypes=[C.c_uint,U8,C.c_uint,U64]
 rows=[]; edges=outputs=0
 for p in range(4,23):
  eps=F(49,132*(1<<p));nodes=get_nodes(p);ids={r['word']:i for i,r in enumerate(nodes)}
  n=len(nodes);depth=max(len(r['word']) for r in nodes)
  nc=C.c_uint64();dc=C.c_uint();assert lib.cf_count(p,C.byref(nc),C.byref(dc))==0
  assert (nc.value,dc.value)==(n,depth)
  B=(n-1).bit_length();BP=depth+1; nb=(B+7)//8; nbp=(BP+7)//8
  states={k:(C.c_uint8*z)() for k,z in [('rank',nb),('path',nbp)]};out=(C.c_uint64*4)()
  for i,r in enumerate(nodes):
   w=r['word'];center=(r['lo']+r['hi'])/2 if r['leaf'] else value(w,F(1,2))
   assert max(abs(r['lo']-center),abs(r['hi']-center))<=2*eps if r['leaf'] else True
   for v,bts,index in [('rank',nb,i),('path',nbp,int('1'+w,2))]:
    state=states[v]
    for k,b in enumerate(index.to_bytes(bts,'little')):state[k]=b
    assert getattr(lib,f'cf_{v}_output')(p,state,bts,out)==0
    numerator=out[0]+(out[1]<<64);denominator=out[2]+(out[3]<<64)
    assert F(numerator,denominator)==F(3,4)-center/2
    assert 0<numerator<denominator;outputs+=1
    for b in (0,1):
     raw=str(b)+w;new=raw
     for k in range(1,len(raw)+1):
      l,h=bounds(raw[:k])
      if h-l<=4*eps:new=raw[:k];break
     assert new in ids
     nr=nodes[ids[new]]
     if r['leaf']:
      assert nr['leaf']
      assert nr['lo']<=bayes(r['lo'],b)<=bayes(r['hi'],b)<=nr['hi']
     else:
      true=bayes(value(w,F(1,2)),b)
      assert nr['lo']<=true<=nr['hi']
      if not nr['leaf']:assert true==value(new,F(1,2))
     for k,bb in enumerate(index.to_bytes(bts,'little')):state[k]=bb
     assert getattr(lib,f'cf_{v}_step')(p,state,bts,b)==0
     expected=ids[new] if v=='rank' else int('1'+new,2)
     assert int.from_bytes(bytes(state),'little')==expected;edges+=1
  previous=json.loads((ROOT/'reference'/f'certificate_p{p:02d}.json').read_text())
  xs=[value(w[::-1],F(1,2)) for w in previous['witness_words']]
  assert all((b-a)/2>2*eps for a,b in zip(xs,xs[1:]))
  LB=(len(xs)-1).bit_length()
  row={'p':p,'epsilon':str(eps),'nodes':n,'leaves':(n+1)//2,'depth':depth,
       'lower_states':len(xs),'lower_bits':LB,'rank_bits':B,'path_bits':BP,
       'rank_bytes':nb,'path_bytes':nbp,'rank_exact':B==LB,'path_exact':BP==LB,
       'old_nodes':previous['upper_states'],'old_table_bytes':16+8*previous['upper_states']}
  rows.append(row);print(json.dumps(row),flush=True)
 # Global finite-width proof for allowed API p<=30.
 assert F(16,49)**20/F(4)<=F(49,33*(1<<30))
 bound=7**20
 assert bound<2**64 and 6*bound<2**127
 assert 49*(16*bound)**2<2**128
 assert 4*(16*bound)**2<2**128
 assert 528*9**20*(1<<30)<2**128
 result={'passed':True,'configs':rows,'C_transition_checks':edges,'C_output_checks':outputs,
         'rank_exact_p':[r['p'] for r in rows if r['rank_exact']],
         'path_exact_p':[r['p'] for r in rows if r['path_exact']],
         'finite_width_proof':{'max_supported_p':30,'max_depth_bound':20,'matrix_entry_bound':bound,
          'largest_upper_bound':49*(16*bound)**2,'integer_type':'uint64 coordinates, unsigned __int128 products'},
         'runtime_tables':0,'immutable_parameter_payload_bytes':2}
 (OUT/'COMPACT_VERIFIED.json').write_text(json.dumps(result,indent=2)+'\n')
 print('DONE',edges,outputs,flush=True)
if __name__=='__main__':main()
