"""Integer specification of the finite C probability oracle; no floating point."""
Q=1<<24
P=1<<16
LN2Q=11629079
def trunc(a,b):
    return a//b if a>=0 else -((-a)//b)
def weight(t,b):
    den=5*Q+2*t if b else 7*Q-2*t
    root=3*Q*Q//den
    child=(3 if b else 1)*Q*(Q+t)//den
    z=trunc((den-6*Q)*Q,den+6*Q)
    z2=z*z//Q;z3=trunc(z2*z,Q);z5=trunc(z3*z2,Q)
    logv=LN2Q+2*(z+trunc(z3,3)+trunc(z5,5))
    x=2*1003*logv//200000
    e=Q-x+x*x//(2*Q)
    a=root*e//Q
    v=trunc(child-Q//2,3);v2=v*v//Q;v4=v2*v2//Q
    H=Q+99577*v2//100000+125691*v4//100000
    return a*H//Q
def branch(S,T):
    t=S*Q//(S+T);w0=weight(t,0);w1=weight(t,1)
    return w0*P//(w0+w1)
def bits(p):return (101003*p+199999)//200000+3
def append(z,b):
    a,c=z[0],z[2]
    return (6*a+3*z[1],a+2*z[1],6*c+3*z[3],c+2*z[3]) if b else (2*a+z[1],3*a+6*z[1],2*c+z[3],3*c+6*z[3])
def leaf(z,p):
    a,b,c,d=z;dl=3*(a+c)+5*(b+d);dr=5*(a+c)+3*(b+d)
    return 528*(1<<p)*(a*d-b*c)<=49*dl*dr
def encode(p,w):
    lo,N,z=0,1<<bits(p),(1,0,0,1)
    for b in map(int,w):
        assert not leaf(z,p)
        k=branch(z[0]+z[2],z[1]+z[3]);left=(N-1)*k//P
        assert 0<left<N-1
        lo+=1+(left if b else 0);N=N-1-left if b else left
        z=append(z,b)
    return lo,N
