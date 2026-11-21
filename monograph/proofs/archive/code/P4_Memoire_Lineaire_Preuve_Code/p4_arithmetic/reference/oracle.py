from fractions import Fraction as F
# Independent projective belief coordinates, derived from original T_b.
T=[((F(1,2),F(1,4)),(F(1,12),F(1,6))),((F(1,6),F(1,12)),(F(1,4),F(1,2)))]
RB=[]
for ((a,b),(c,d)) in T:
    vals=[12*(d-b),12*b,12*(c+d-a-b),12*(a+b)]
    assert all(x.denominator==1 for x in vals)
    RB.append(tuple(int(x) for x in vals))
assert RB==[(-1,3,-6,9),(5,1,6,3)]
def bayes(x,b):
    ((a,c),(d,e))=T[b]
    left=(1-x)*a+x*d;right=(1-x)*c+x*e
    return right/(left+right)
def nodes(p):
    eps=F(49,132*(1<<p));out=[]
    def visit(w,a,b,c,d):
        lo=F(3*a+8*b,3*c+8*d);hi=F(5*a+8*b,5*c+8*d)
        leaf=hi-lo<=4*eps
        center=(lo+hi)/2 if leaf else F(a+2*b,c+2*d)
        out.append({'word':w,'lo':lo,'hi':hi,'leaf':leaf,'center':center,'out':F(3,4)-center/2})
        if not leaf:
            for bit,(A,B,C,D) in enumerate(RB):
                visit(w+str(bit),a*A+b*C,a*B+b*D,c*A+d*C,c*B+d*D)
    visit('',1,0,0,1);return out
