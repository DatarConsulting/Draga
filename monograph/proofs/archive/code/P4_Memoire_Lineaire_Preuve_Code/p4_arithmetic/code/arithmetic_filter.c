#include <stdint.h>
#include <stddef.h>
typedef unsigned __int128 u128;
typedef struct { uint64_t a,b,c,d; unsigned depth; uint64_t path; } Cursor;
/* Fixed source only. p and bytes are immutable caller parameters. No table,
   allocation, mutable global, or external clock. p in [4,30], depth <=20; exact rational output. */
static Cursor root(void) { Cursor z={1,0,0,1,0,1}; return z; }
static Cursor append(Cursor z,unsigned b) {
    uint64_t a=z.a,c=z.c;
    if(b==0) { z.a=2*a+z.b; z.b=3*a+6*z.b; z.c=2*c+z.d; z.d=3*c+6*z.d; }
    else { z.a=6*a+3*z.b; z.b=a+2*z.b; z.c=6*c+3*z.d; z.d=c+2*z.d; }
    z.path=(z.path<<1)|b; z.depth++; return z;
}
static int leaf(Cursor z,unsigned p) {
    u128 det=(u128)z.a*z.d-(u128)z.b*z.c;
    u128 dl=3*((u128)z.a+z.c)+5*((u128)z.b+z.d);
    u128 dr=5*((u128)z.a+z.c)+3*((u128)z.b+z.d);
    return (det*528*((uint64_t)1<<p))<=49*dl*dr;
}
static uint64_t load(const uint8_t *s,unsigned nb) {
    uint64_t v=0;for(unsigned i=0;i<nb;i++)v|=(uint64_t)s[i]<<(8*i);return v;
}
static void save(uint8_t *s,unsigned nb,uint64_t v) {for(unsigned i=0;i<nb;i++)s[i]=(uint8_t)(v>>(8*i));}
/* Binary probability P(0) returned as 128-bit numerator and denominator,
   each split into low/high 64-bit limbs. P(1)=1-P(0). */
static void output(Cursor z,unsigned p,uint64_t out[4]) {
    u128 num,den;
    if(leaf(z,p)) {
        u128 nl=3*(u128)z.a+5*(u128)z.b,nr=5*(u128)z.a+3*(u128)z.b;
        u128 dl=nl+3*(u128)z.c+5*(u128)z.d,dr=nr+5*(u128)z.c+3*(u128)z.d;
        den=4*dl*dr;num=3*dl*dr-nl*dr-nr*dl;
    } else {
        u128 n=(u128)z.a+z.b,S=n+z.c+z.d;den=4*S;num=3*S-2*n;
    }
    out[0]=(uint64_t)num;out[1]=(uint64_t)(num>>64);out[2]=(uint64_t)den;out[3]=(uint64_t)(den>>64);
}
#define Q ((int64_t)1<<24)
#define P ((uint64_t)1<<16)
static uint64_t weight(int64_t t,unsigned b) {
    int64_t den=b?5*Q+2*t:7*Q-2*t;
    int64_t root=3*Q*Q/den,child=(b?3:1)*Q*(Q+t)/den;
    int64_t z=(den-6*Q)*Q/(den+6*Q),z2=z*z/Q,z3=z2*z/Q,z5=z3*z2/Q;
    int64_t logv=11629079+2*(z+z3/3+z5/5);
    int64_t x=2*1003*logv/200000,e=Q-x+x*x/(2*Q),a=root*e/Q;
    int64_t v=(child-Q/2)/3,v2=v*v/Q,v4=v2*v2/Q;
    int64_t H=Q+99577*v2/100000+125691*v4/100000;
    return (uint64_t)(a*H/Q);
}
uint32_t ac_branch(uint64_t S,uint64_t T) {
    /* Exact floor(2^24*S/(S+T)), without a 128-bit division helper. */
    uint64_t den=S+T,rem=S,t=0;
    for(unsigned i=0;i<24;i++){rem*=2;t<<=1;if(rem>=den){rem-=den;t++;}}
    uint64_t w0=weight((int64_t)t,0),w1=weight((int64_t)t,1);
    return (uint32_t)(w0*P/(w0+w1));
}
unsigned ac_bits(unsigned p){return (101003*p+199999)/200000+3;}
static int split(Cursor z,uint64_t *lo,uint64_t *N,unsigned b){
    uint64_t left=(*N-1)*ac_branch(z.a+z.c,z.b+z.d)/P;
    if(!left||left>=*N-1)return -1;
    *lo+=1+(b?left:0);*N=b?*N-1-left:left;return 0;
}
int ac_encode(unsigned p,uint64_t path,uint64_t *id){
    if(p<4||p>30||!path)return -1;
    unsigned depth=0;for(uint64_t v=path;v>1;v>>=1)depth++;
    if(depth>20)return -1;
    Cursor z=root();uint64_t lo=0,N=(uint64_t)1<<ac_bits(p);
    for(unsigned i=depth;i>0;i--){
        unsigned b=(path>>(i-1))&1;
        if(leaf(z,p)||split(z,&lo,&N,b))return -2;
        z=append(z,b);
    }
    *id=lo;return 0;
}
static int decode(unsigned p,uint64_t id,Cursor *out){
    if(p<4||p>30)return -1;
    Cursor z=root();uint64_t lo=0,N=(uint64_t)1<<ac_bits(p);
    if(id>=N)return -1;
    while(id!=lo){
        if(leaf(z,p)||z.depth>=20)return -2;
        uint64_t left=(N-1)*ac_branch(z.a+z.c,z.b+z.d)/P;
        if(!left||left>=N-1)return -3;
        unsigned b=id>=lo+1+left;
        lo+=1+(b?left:0);N=b?N-1-left:left;z=append(z,b);
    }
    *out=z;return 0;
}
int ac_decode(unsigned p,uint64_t id,uint64_t *path){Cursor z;int rc=decode(p,id,&z);if(!rc)*path=z.path;return rc;}
int ac_output(unsigned p,const uint8_t *state,unsigned nb,uint64_t out[4]){
    if(p<4||p>30||nb<1||nb>4||8*nb<ac_bits(p))return -1;
    Cursor z;if(decode(p,load(state,nb),&z))return -2;output(z,p,out);return 0;
}
int ac_advance(unsigned p,uint8_t *state,unsigned nb,unsigned b,uint64_t out[4]){
    if(p<4||p>30||nb<1||nb>4||8*nb<ac_bits(p)||b>1)return -1;
    Cursor old;if(decode(p,load(state,nb),&old))return -2;
    Cursor z=root();uint64_t lo=0,N=(uint64_t)1<<ac_bits(p);
    if(split(z,&lo,&N,b))return -3;z=append(z,b);
    for(unsigned i=old.depth;i>0&&!leaf(z,p);i--){
        b=(old.path>>(i-1))&1;if(split(z,&lo,&N,b))return -3;z=append(z,b);
    }
    save(state,nb,lo);output(z,p,out);return 0;
}
