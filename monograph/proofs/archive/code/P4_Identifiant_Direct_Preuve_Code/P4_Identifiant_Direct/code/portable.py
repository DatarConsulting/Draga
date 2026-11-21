"""Exact arbitrary-precision direct suffix filter.

Persistent state is ONE integer, initially 1. p is an immutable parameter.
For all integer p>=4: polynomial bit time, O(p) bits of working/persistent
memory. Optimal bit budgets are only claimed on the certified finite grid.
No tree enumeration, recursion, lookup table, or hidden history.
"""
from fractions import Fraction

def fold(p,path):
    if p not in (19,21):return path
    B=(p+1)//2;H=1<<(B-1);v=(H-1)//3
    if path>>1==H+v:return H+(path&1)
    if path>>1==2*H-1-v:return H+2+(path&1)
    return path

def unfold(p,index):
    if p not in (19,21):return index
    B=(p+1)//2;H=1<<(B-1);v=(H-1)//3
    if not H<=index<H+4:return index
    offset=index-H
    parent=H+(H-1-v if offset>>1 else v)
    return (parent<<1)|(offset&1)

def advance(p,index,bit):
    """Precondition: valid reachable index, p>=4, bit in {0,1}.

    Returns (new_index, exact_P0, diagnostics). Diagnostics are transient.
    """
    if p<4 or bit not in (0,1) or index<1:raise ValueError('invalid argument')
    path=unfold(p,index);depth=path.bit_length()-1
    # Construct newest bit followed by the known suffix, excluding sentinel.
    raw=(bit<<depth)|(path^(1<<depth))
    A,B,C,D=1,0,0,1;newpath=1;processed=0
    for i in range(depth,-1,-1):
        b=(raw>>i)&1
        if b==0:A,B,C,D=2*A+B,3*A+6*B,2*C+D,3*C+6*D
        else:A,B,C,D=6*A+3*B,A+2*B,6*C+3*D,C+2*D
        newpath=(newpath<<1)|b;processed+=1
        determinant=A*D-B*C
        dl=3*(A+C)+5*(B+D);dr=5*(A+C)+3*(B+D)
        isleaf=528*(1<<p)*determinant<=49*dl*dr
        if isleaf:break
    if isleaf:
        nl=3*A+5*B;nr=5*A+3*B
        numerator=3*dl*dr-nl*dr-nr*dl;denominator=4*dl*dr
    else:
        s=A+B+C+D;n=A+B
        numerator=3*s-2*n;denominator=4*s
    return fold(p,newpath),Fraction(numerator,denominator),{
        'processed_bits':processed,'coordinate_bits':max(A,B,C,D).bit_length(),
        'output_bits':max(numerator,denominator).bit_length()}
