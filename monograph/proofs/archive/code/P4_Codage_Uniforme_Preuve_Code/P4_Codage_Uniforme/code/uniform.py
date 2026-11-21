"""Uniform, exact integer construction. Immutable DP data is explicitly charged.
The Python backend uses arbitrary-size integers. The C backend is p=4..30.
No suffix tree is enumerated during construction, encoding or decoding.
"""
from pathlib import Path
from array import array
import ctypes as C
import subprocess
from fractions import Fraction

ROOT=Path(__file__).resolve().parents[1]
def exact_update(p,word,bit):
    """Update the actual stopped suffix, not the quantized coding state."""
    if bit not in (0,1):raise ValueError('binary symbol required')
    a,b,c,d=1,0,0,1;new=[]
    for char in str(bit)+word:
        if char=='0':a,b,c,d=2*a+b,3*a+6*b,2*c+d,3*c+6*d
        else:a,b,c,d=6*a+3*b,a+2*b,6*c+3*d,c+2*d
        new.append(char)
        dl=3*(a+c)+5*(b+d);dr=5*(a+c)+3*(b+d)
        leaf=528*(1<<p)*(a*d-b*c)<=49*dl*dr
        if leaf:break
    if leaf:
        nl=3*a+5*b;nr=5*a+3*b
        out=Fraction(3*dl*dr-nl*dr-nr*dl,4*dl*dr)
    else:out=Fraction(3*(a+b+c+d)-2*(a+b),4*(a+b+c+d))
    return ''.join(new),out
class Context(C.Structure):
    _fields_=[('p',C.c_uint),('m',C.c_uint),('J',C.c_uint),('R',C.c_uint),
              ('next',C.POINTER(C.c_uint32)),('cost',C.POINTER(C.c_uint32)),('counts',C.POINTER(C.c_uint64))]

def compile_lib():
    dest=ROOT/'code/libuniform.so'
    subprocess.run(['gcc','-std=c11','-O2','-Wall','-Wextra','-Wno-misleading-indentation',
                    '-shared','-fPIC',str(ROOT/'code/uniform_filter.c'),'-o',str(dest)],check=True)
    return dest

def ceil_log_ratio(n,d):
    assert n>d>0
    k=max(0,n.bit_length()-d.bit_length())
    if d<<k<n:k+=1
    assert (d<<k)>=n and (k==0 or (d<<(k-1))<n)
    return k

def weights(p):
    if p<4:raise ValueError('p >= 4 required')
    m=8*(p+1);J=3*m//4+1;R=m*p+(m+5)//6
    nxt=[];cost=[];den=(3*m)**(2*m)
    for j in range(m//8,7*m//8+1):
        for b in (0,1):
            D=7*m-2*j if b==0 else 5*m+2*j
            numerator=m*(m+j)*(1 if b==0 else 3)
            dest=(2*numerator+D)//(2*D)
            assert m//8<=dest<=7*m//8
            nxt.append(dest-m//8)
            k=ceil_log_ratio(D**(2*m),den)
            assert m<=k<=3*m
            cost.append(k)
    return m,J,R,nxt,cost

class Codec:
    def __init__(self,p,backend='python',lib=None):
        self.p=p;self.m,self.J,self.R,nxt,cost=weights(p)
        self.next=nxt;self.cost=cost;self.initial=3*self.m//8
        J,R=self.J,self.R
        self.backend=backend
        if backend=='c':
            if not 4<=p<=30:raise ValueError('C backend supports p=4..30')
            self.lib=C.CDLL(str(lib or compile_lib()))
            self.lib.uc_build.argtypes=[C.c_uint,C.c_uint,C.c_uint,C.POINTER(C.c_uint32),C.POINTER(C.c_uint32),C.POINTER(C.c_uint64)]
            self.lib.uc_encode.argtypes=[C.POINTER(Context),C.c_uint64,C.POINTER(C.c_uint64)]
            self.lib.uc_decode.argtypes=[C.POINTER(Context),C.c_uint64,C.POINTER(C.c_uint64)]
            self.lib.uc_advance.argtypes=[C.POINTER(Context),C.POINTER(C.c_uint8),C.c_uint,C.c_uint,C.POINTER(C.c_uint64)]
            self.lib.uc_output.argtypes=[C.POINTER(Context),C.POINTER(C.c_uint8),C.c_uint,C.POINTER(C.c_uint64)]
            self.na=(C.c_uint32*len(nxt))(*nxt);self.ka=(C.c_uint32*len(cost))(*cost)
            self.counts=array('Q',[0])*(J*(R+1))
            self.ca=(C.c_uint64*len(self.counts)).from_buffer(self.counts)
            rc=self.lib.uc_build(p,J,R,self.na,self.ka,self.ca)
            if rc:raise RuntimeError(('C build failed',rc))
            self.ctx=Context(p,self.m,J,R,self.na,self.ka,self.ca)
        elif backend=='python':
            self.counts=[0]*(J*(R+1))
            for r in range(R+1):
                for j in range(J):
                    n=1
                    for b in (0,1):
                        k=cost[2*j+b]
                        if r>=k:n+=self.counts[(r-k)*J+nxt[2*j+b]]
                    self.counts[r*J+j]=n
        else:raise ValueError('unknown backend')
        self.size=self.counts[R*J+self.initial];self.bits=(self.size-1).bit_length()
        self.bytes=(self.bits+7)//8
    def encode(self,word):
        j=self.initial;r=self.R;rank=0
        for char in word:
            b=int(char);k=self.cost[2*j+b]
            if k>r:raise ValueError('outside the coding language')
            rank+=1
            k0=self.cost[2*j]
            if b and r>=k0:rank+=self.counts[(r-k0)*self.J+self.next[2*j]]
            r-=k;j=self.next[2*j+b]
        return rank
    def decode(self,index):
        if not 0<=index<self.size:raise ValueError('invalid ID')
        j=self.initial;r=self.R;word=[]
        while index:
            index-=1;k0=self.cost[2*j]
            left=0 if r<k0 else self.counts[(r-k0)*self.J+self.next[2*j]]
            b=int(index>=left)
            if b:index-=left
            k=self.cost[2*j+b]
            assert k<=r
            r-=k;j=self.next[2*j+b];word.append(str(b))
        return ''.join(word)
    def word_cost(self,word):
        j=self.initial;s=0
        for char in word:
            b=int(char);s+=self.cost[2*j+b];j=self.next[2*j+b]
        return s
    def advance(self,index,bit):
        word,out=exact_update(self.p,self.decode(index),bit)
        return self.encode(word),out

def count_rolling(p,start,budget):
    """Exact count with only max(cost)+1 budget rows; no suffix enumeration.
    Can replace full table queries at a polynomial extra time cost.
    """
    m,J,R,nxt,cost=weights(p);H=max(cost)+1
    rows=[[0]*J for _ in range(H)]
    for r in range(budget+1):
        row=rows[r%H]
        for j in range(J):
            n=1
            for b in (0,1):
                k=cost[2*j+b]
                if r>=k:n+=rows[(r-k)%H][nxt[2*j+b]]
            row[j]=n
    return rows[budget%H][start]

class RollingCodec:
    """Same IDs, no persistent full count table; each count is recomputed.
    Fixed grid data O(p log p) bits, transient counters O(p^3) bits.
    """
    def __init__(self,p):
        self.p=p;self.m,self.J,self.R,self.next,self.cost=weights(p)
        self.initial=3*self.m//8
        self.size=self.count(self.initial,self.R);self.bits=(self.size-1).bit_length()
    def count(self,start,budget):
        J=self.J;H=max(self.cost)+1
        rows=[[0]*J for _ in range(H)]
        for r in range(budget+1):
            row=rows[r%H]
            for j in range(J):
                n=1
                for b in (0,1):
                    k=self.cost[2*j+b]
                    if r>=k:n+=rows[(r-k)%H][self.next[2*j+b]]
                row[j]=n
        return rows[budget%H][start]
    def encode(self,word):
        j=self.initial;r=self.R;rank=0
        for char in word:
            b=int(char);k=self.cost[2*j+b];k0=self.cost[2*j]
            if k>r:raise ValueError('outside language')
            rank+=1
            if b and r>=k0:rank+=self.count(self.next[2*j],r-k0)
            r-=k;j=self.next[2*j+b]
        return rank
    def decode(self,index):
        if not 0<=index<self.size:raise ValueError('invalid ID')
        j=self.initial;r=self.R;word=[]
        while index:
            index-=1;k0=self.cost[2*j]
            left=0 if r<k0 else self.count(self.next[2*j],r-k0)
            b=int(index>=left)
            if b:index-=left
            k=self.cost[2*j+b];assert k<=r
            r-=k;j=self.next[2*j+b];word.append(str(b))
        return ''.join(word)
    def advance(self,index,bit):
        word,out=exact_update(self.p,self.decode(index),bit)
        return self.encode(word),out
