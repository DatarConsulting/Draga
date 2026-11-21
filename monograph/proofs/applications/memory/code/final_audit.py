"""Cross-check deployed packets, target rank, and the local-smoothing KL bound."""
import json,struct,hashlib,math
from pathlib import Path
from fractions import Fraction as F
import numpy as np
import verify as v
ROOT=Path(__file__).resolve().parents[1];D=ROOT/'results'
def upper_log(x):
    count=0
    while x>=2:x=x/2;count+=1
    def atanh_upper(t):
        u=(t-1)/(t+1);power=u;total=F(0)
        for j in range(24):total+=2*power/F(2*j+1);power*=u*u
        return total+2*power/(49*(1-u*u))
    return count*atanh_upper(F(2))+atanh_upper(x)
def run():
    packets=0;targets=0
    for p in sorted((ROOT/'models').glob('*.json')):
        m=json.loads(p.read_text());A=[v.mat(a) for a in m['reduced_A']];al=v.mat(m['reduced_alpha']);be=v.mat(m['reduced_beta'])
        assert any(v.det(np.concatenate([al,al@a],axis=0)) for a in A)
        assert any(v.det(np.concatenate([be,a@be],axis=1)) for a in A)
        assert not np.array_equal(A[0]@A[1],A[1]@A[0]);targets+=1
        for policy in ['uniform','mixed']:
            c=m['best'][policy];f=m['frames'][c['basis']];raw=(ROOT/'models'/f"{m['id']}_{policy}.bin").read_bytes()
            assert raw[:6]==b'R2Q1\x02\x18' and list(raw[6:8])==c['bits']
            assert list(struct.unpack('<2i',raw[8:16]))==c['shifts']
            ai=[int(F(x)*2**24) for a in f['A'] for row in a for x in row]
            assert list(struct.unpack('<8i',raw[16:48]))==ai
            assert list(struct.unpack('<2i',raw[48:56]))==[int(F(x)*2**24) for row in f['alpha'] for x in row]
            assert list(struct.unpack('<2i',raw[56:64]))==[int(F(x[0])*2**24) for x in f['beta']]
            packed=int.from_bytes(raw[64:],'little');off=0
            for pi,x in zip(c['bits'],c['initial_integers']):
                val=(packed>>off)&(2**pi-1);assert val==x%2**pi;off+=pi
            rec=json.loads((D/f"replay_{m['id']}_{policy}.json").read_text());assert hashlib.sha256(raw).hexdigest()==rec['packet_sha256'];packets+=1
    prob=json.loads((D/'probability.json').read_text());models={c['id']:c for c in prob['synthetic']+prob['trained']}
    kl=json.loads((D/'KL_realization.json').read_text());checked=0
    for row in kl['records']:
        c=models[row['id']];d=F(row['d_for_KL']);lam=F(row['lambda']);tau=F(row['tail_stopping']);m=c['alphabet_size'];N=row['tail_cutoff']
        assert d>=F(c['total_word_TV_bound']) and 0<d<1 and 0<lam<1
        assert F(row['smoothed_TV_bound'])==F(c['total_word_TV_bound'])+lam/F(c['stopping_min'])<=F(1,100)
        psi=d*N+(1-tau)**(N+1)/tau;assert psi==F(row['tail_first_moment_bound'])
        U=upper_log(1/(1-d))+upper_log(1/(1-lam))/tau+d*upper_log(2/lam)+psi*upper_log(2*m/lam)
        assert F(row['KL_rational_upper'])>=U
        assert row['state_bytes']==c['state_bytes'] and row['maximum_unsigned_arithmetic_bits']<=64
        checked+=1
    report={'all_passed':True,'rank_two_targets_and_noncommutativity':targets,'binary_packets_verified':packets,'KL_certificates_verified':checked,
      'note':'Independent packet decoding and rational re-evaluation; compiler and theorem proof remain ordinary mathematical review, not formal proof-assistant verification.'}
    (D/'final_audit.json').write_text(json.dumps(report,indent=2));print(report)
if __name__=='__main__':run()
