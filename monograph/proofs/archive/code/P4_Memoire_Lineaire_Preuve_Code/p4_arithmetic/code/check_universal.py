import json,sys,time
from pathlib import Path
from fractions import Fraction as F
from universal import Codec,Dyadic,horizon,transfer
import universal
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'reference'))
from oracle import nodes
def main():
    start=time.perf_counter();rows=[]
    for p in [4,6]:
        metrics={'attempts':0,'max_precision':0}
        def tracked(q):
            metrics['attempts']+=1;metrics['max_precision']=max(q,metrics['max_precision'])
            return Dyadic(q)
        # Instrumentation belongs solely to this verifier, never to the filter.
        universal.Dyadic=tracked
        c=Codec(p);ns=nodes(p);byword={v['word']:v for v in ns};ids={w:c.encode(w)[0] for w in byword}
        assert len(set(ids.values()))==len(ns)
        edges=0
        for w,i in ids.items():
            assert c.decode(i)==w
            assert c.output(i)==byword[w]['out']
            for b in [0,1]:
                raw=str(b)+w;new=raw
                for k in range(1,len(raw)+1):
                    if byword[raw[:k]]['leaf']:new=raw[:k];break
                assert c.advance(i,b)==ids[new];edges+=1
        rows.append({'p':p,'states':len(ns),'bits':c.B,'horizon':c.n,'transitions':edges,'exact_outputs':len(ns),
                     'oracle_interval_attempts':metrics['attempts'],
                     'summed_terms':metrics['attempts']*2*(1<<c.n),'max_precision':metrics['max_precision']})
        print(rows[-1],flush=True)
    universal.Dyadic=Dyadic
    # Independent high precision arithmetic checks of 64-bit interval enclosure,
    # not used as proof; the enclosure proof comes from the explicit series tails.
    import mpmath as mp
    mp.mp.dps=90
    c=Dyadic(64);checks=0
    for a,b in [(1,1),(2,1),(81,16),(10**60+1,7**30),(7**100,9**30)]:
        lo,hi=c.log_ratio(a,b);v=mp.log(mp.mpf(a)/b)
        assert mp.mpf(lo)/c.S<=v<=mp.mpf(hi)/c.S;checks+=1
    # Compute a nontrivial branch using a very long exact projective argument.
    # p remains 4 to isolate input integer length from transfer horizon.
    cc=Codec(4);long=cc.branch(10**60+7,10**60-3)
    assert long==cc.Q//2
    out={'passed':True,'graphs':rows,'independent_log_checks':checks,
         'long_integer_branch_checked_bits':(10**60+7).bit_length(),
         'large_p_root_refinement_executed':False,'seconds_local':time.perf_counter()-start}
    (ROOT/'results/UNIVERSAL.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
