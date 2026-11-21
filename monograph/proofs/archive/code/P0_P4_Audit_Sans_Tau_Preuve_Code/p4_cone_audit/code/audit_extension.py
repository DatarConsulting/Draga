"""Executed audit: exact references are confined to this verifier."""
from fractions import Fraction as F
from pathlib import Path
from math import comb
import json,random,time,sys,importlib.util
from certificate import certify,hmm,CertificateError,ceil_log2,ceil
from cases import CASES
from model import Model
from codec import Codec
from dyadic import Dyadic
from oracle import word_center,power,transfer_coefficients,ratio,precision_budget
from verify import run_case,mul,belief,value,output

ROOT=Path(__file__).resolve().parents[1]
EXAMPLES=[('beyond_symmetric',hmm('1/4','1/4','6/25','19/25'),(F(1,4),F(3,4))),
          ('beyond_asymmetric',hmm('6/25','13/50','6/25','3/4'),(F(1,4),F(3,4)))]

def save(name,data):
    (ROOT/'results'/name).write_text(json.dumps(data,indent=2)+'\n')

def geometry_audit():
    old_path=ROOT.parent/'p4_parametric/code/certificate.py'
    # Bundled copy used when the independent reference directory is absent.
    if not old_path.exists():old_path=ROOT/'reference/old_certificate.py'
    spec=importlib.util.spec_from_file_location('old_certificate',old_path)
    old=importlib.util.module_from_spec(spec);sys.modules[spec.name]=old;spec.loader.exec_module(old)
    rows=[]
    for name,ms,J in EXAMPLES:
        c=certify(ms,J=J)
        try:old.certify(ms,J=J);raise AssertionError('old certificate must reject')
        except old.CertificateError as e:assert e.condition=='transfer_contraction'
        assert c.tau>1 and 0<c.cone_contraction<1 and c.cone_a>c.cone_b
        assert c.cone_contraction**c.cone_half_life<=F(1,2)
        rows.append({'name':name,'old_tau':str(c.tau),'old_tau_decimal':float(c.tau),
                     'new_gamma_decimal':float(c.cone_contraction),'half_life':c.cone_half_life,
                     'rho':float(c.rho),'dual_rho':float(c.grho),'gap':float(c.separation),
                     'coarse_time_degree':c.cone_half_life+5,'certificate':c.json()})
    save('NEW_CLASS.json',{'passed':True,'models':rows})
    return rows

def graph_audit():
    rows=[]
    schedule=[(case,p) for case in EXAMPLES for p in (4,8,12)]+[(CASES[1],8)]
    for case,p in schedule:
        row=run_case(case,p);rows.append(row)
        print(json.dumps({k:row[k] for k in ('name','p','states','transitions','bits','seconds')}),flush=True)
        save('GRAPHS.json',{'passed':len(rows)==len(schedule),'configs':rows,
             'total_states':sum(x['states'] for x in rows),'total_transitions':sum(x['transitions'] for x in rows),
             'total_outputs':sum(x['outputs'] for x in rows),'max_TV_over_epsilon':max(x['max_TV_over_epsilon_decimal'] for x in rows)})
    return rows

def exact_jet_audit():
    rng=random.Random(503172);word_checks=0;coefficient_checks=0;power_checks=0
    for case in EXAMPLES+[CASES[1]]:
        cert=certify(case[1],J=case[2]);model=Model.from_rationals(16,cert)
        rad=max(F(1,2)-cert.Jdual[0],cert.Jdual[1]-F(1,2))
        for n in (0,1,2,4,8,16,32):
            for _ in range(8):
                w=rng.randrange(1<<n);c=Dyadic(256);D,u=word_center(model,w,n,c)
                P=(1,0,0,1)
                for j in range(n-1,-1,-1):P=mul(P,model.mats[(w>>j)&1])
                a,b,cc,d=P;den=F(a+b+cc+d,2)
                exactD=F(a*d-b*cc,den**2)
                exactu=F(a+b-cc-d,den)
                assert F(D[0],c.S)<=exactD<=F(D[1],c.S)
                assert F(u[0],c.S)<=exactu<=F(u[1],c.S)
                assert abs(exactu*rad)<=cert.taylor_radius
                word_checks+=1
        # sigma=1 yields rational powers and exact Taylor coefficients.
        for n in range(5):
            degree=12;c=Dyadic(160)
            coeffs=transfer_coefficients(model,n,F(1),degree,c)
            refs=[F(0)]*(degree+1)
            for w in range(1<<n):
                P=(1,0,0,1)
                for j in range(n-1,-1,-1):P=mul(P,model.mats[(w>>j)&1])
                a,b,cc,d=P;den=F(a+b+cc+d,2)
                exactD=F(a*d-b*cc,den**2);v=-F(a+b-cc-d,den)*rad
                for k in range(degree+1):refs[k]+=exactD*(k+1)*v**k
            for enclosure,ref in zip(coeffs,refs):
                assert F(enclosure[0],c.S)<=ref<=F(enclosure[1],c.S)
                coefficient_checks+=1
    # Exact integer comparisons verify dyadic fractional powers without floats.
    for q in (24,48,96):
        c=Dyadic(q)
        for den in (1,2,4,8,16,32):
            for numerator in (1,3,7,19):
                s=F(numerator,den)
                for value_ in (F(1,1024),F(1,3),F(7,8),F(1)):
                    z=power(c,c.rat(value_),s)
                    assert F(z[0],c.S)**s.denominator<=value_**s.numerator<=F(z[1],c.S)**s.denominator
                    power_checks+=1
    data={'passed':True,'word_D_and_u_checks':word_checks,'exact_Taylor_coefficients':coefficient_checks,
          'exact_power_checks':power_checks,'scope':'All comparisons rational; verifier uses exact long products only.'}
    save('EXACT_ORACLE_AUDIT.json',data);print(json.dumps(data),flush=True)

def numeric_global_check():
    import mpmath as mp
    from check_oracle_reference import V,real
    mp.mp.dps=90;checks=0;largest=0.0
    configs=json.loads((ROOT/'results/GRAPHS.json').read_text())['configs']
    for row in configs:
        case=next(x for x in EXAMPLES+CASES if x[0]==row['name'])
        H=certify(case[1],J=case[2]);model=Model.from_rationals(row['p'],H)
        n=row['oracle']['horizon'];s=F(row['oracle']['sigma']);bound=F(row['oracle']['max_global_ratio_deviation'])
        for j in range(17):
            t=H.Jdual[0]+F(j,16)*(H.Jdual[1]-H.Jdual[0])
            ref=V(model,real(t),n+1,real(s))/V(model,real(t),n,real(s))
            assert abs(ref-1)<=real(bound)
            largest=max(largest,float(abs(ref-1)/real(bound)));checks+=1
    save('NUMERIC_GLOBAL_REFERENCE.json',{'passed':True,'checks':checks,'decimal_digits':90,
          'max_observed_fraction_of_global_bound':largest,'formal_proof':False})

def stream_audit():
    rows=[]
    for case,p in [(EXAMPLES[0],16),(EXAMPLES[1],24)]:
        H=certify(case[1],J=case[2]);model=Model.from_rationals(p,H);t0=time.perf_counter();codec=Codec(model)
        prep=time.perf_counter()-t0;state=0;x=F(1,2);maxerr=F(0);rng=random.Random(1700+p)
        statebits=[];t0=time.perf_counter()
        for j in range(128):
            b=(j&1) if j<32 else rng.randrange(2)
            x=value(belief(case[1][b]),x);state=codec.advance(state,b)
            err=abs(F(codec.output(state),model.Q)-output(case[1],x));assert err<=F(1,1<<p)
            maxerr=max(maxerr,err);statebits.append(state.bit_length())
        elapsed=time.perf_counter()-t0
        rows.append({'name':case[0],'p':p,'steps':128,'B':codec.B,'k':model.k,'q':model.q,'oracle_horizon':codec.n,
                     'oracle_precision':codec.certificate['precision_bits'],'Taylor_degree':codec.certificate['taylor_degree'],
                     'max_TV_over_epsilon':float(maxerr*(1<<p)),'preparation_seconds':prep,
                     'seconds_update_and_output_total':elapsed,'milliseconds_per_step':1000*elapsed/128,
                     'timing_scope':'Python, includes verifier exact posterior calculation; no performance claim'})
        print(json.dumps(rows[-1]),flush=True)
    save('STREAMS.json',{'passed':True,'rows':rows})

def main():
    start=time.perf_counter();geometry_audit();exact_jet_audit();graph_audit();numeric_global_check();stream_audit()
    save('RUN_COMPLETE.json',{'passed':True,'elapsed_seconds':time.perf_counter()-start})
if __name__=='__main__':main()
