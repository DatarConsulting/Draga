"""Autonomous arithmetic suffix code. Persistent history is the ID alone."""
from oracle import branch,compile_oracle

class Codec:
    def __init__(self,model,progress=None):
        self.model=model
        self.sigma,self.n,cert=compile_oracle(model,progress)
        self.B=-((-self.sigma.numerator*model.p)//self.sigma.denominator)+8
        # Small fixed certificate summaries are returned for audit separately;
        # no grid, branch table, history counters, or cached word is retained.
        self.certificate=cert
    def probability(self,z):
        return branch(self.model,self.model.dual(z),self.model.Q,self.n,self.sigma)
    def encode(self,w,truncate=False):
        model=self.model;z=model.root();lo,N=0,1<<self.B
        for ch in w:
            if model.leaf(z):
                if truncate:break
                raise ValueError('path extends beyond a leaf')
            b=int(ch);k,Q=self.probability(z);left=(N-1)*k//Q
            assert 0<left<N-1
            lo+=1+(left if b else 0);N=N-1-left if b else left;z=model.step(z,b)
        return lo,N,z
    def decode(self,i):
        if not 0<=i<1<<self.B:raise ValueError('invalid ID')
        model=self.model;z=model.root();lo,N=0,1<<self.B;w=''
        while i!=lo:
            if model.leaf(z):raise ValueError('unused ID')
            k,Q=self.probability(z);left=(N-1)*k//Q
            assert 0<left<N-1
            b=int(i>=lo+1+left);lo+=1+(left if b else 0);N=N-1-left if b else left
            z=model.step(z,b);w+=str(b)
        return w,z
    def advance(self,i,b):
        assert b in (0,1)
        w,_=self.decode(i);return self.encode(str(b)+w,truncate=True)[0]
    def output(self,i):return self.model.output(self.decode(i)[1])
