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
static Cursor update(Cursor old,unsigned b,unsigned p) {
    Cursor z=append(root(),b);
    for(unsigned i=old.depth;i>0&&!leaf(z,p);i--) z=append(z,(old.path>>(i-1))&1);
    return z;
}
static uint64_t load(const uint8_t *s,unsigned nb) {
    uint64_t v=0;for(unsigned i=0;i<nb;i++)v|=(uint64_t)s[i]<<(8*i);return v;
}
static void save(uint8_t *s,unsigned nb,uint64_t v) {for(unsigned i=0;i<nb;i++)s[i]=(uint8_t)(v>>(8*i));}
static int from_path(unsigned p,uint64_t path,Cursor *result) {
    if(path==0)return -1;unsigned depth=0;for(uint64_t t=path;t>1;t>>=1)depth++;
    if(depth>20)return -1;Cursor z=root();
    for(unsigned i=depth;i>0;i--) {if(leaf(z,p))return -1;z=append(z,(path>>(i-1))&1);}
    *result=z;return 0;
}
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
typedef struct {
    unsigned p,m,J,R;
    const uint32_t *next,*cost;
    const uint64_t *counts;
} UC;
/* DP over grid state and integer cost budget; no suffix enumeration. */
int uc_build(unsigned p,unsigned J,unsigned R,const uint32_t *next,
             const uint32_t *cost,uint64_t *counts) {
    if(p<4||p>30||!J)return -1;
    for(unsigned r=0;r<=R;r++)for(unsigned j=0;j<J;j++) {
        uint64_t n=1;
        for(unsigned b=0;b<2;b++) {
            unsigned k=cost[2*j+b],dest=next[2*j+b];
            if(!k||dest>=J)return -2;
            if(r>=k) {
                uint64_t v=counts[(size_t)(r-k)*J+dest];
                if(UINT64_MAX-n<v)return -3;n+=v;
            }
        }
        counts[(size_t)r*J+j]=n;
    }
    return 0;
}
static uint64_t count0(const UC *c,unsigned j,unsigned r) {
    unsigned k=c->cost[2*j];
    return r<k?0:c->counts[(size_t)(r-k)*c->J+c->next[2*j]];
}
int uc_encode(const UC *c,uint64_t path,uint64_t *id) {
    if(!path)return -1;
    unsigned n=0;for(uint64_t v=path;v>1;v>>=1)n++;
    unsigned j=3*c->m/8,r=c->R;uint64_t rank=0;
    for(unsigned i=n;i>0;i--) {
        unsigned b=(unsigned)((path>>(i-1))&1),k=c->cost[2*j+b];
        if(k>r)return -2;
        rank++;if(b)rank+=count0(c,j,r);
        r-=k;j=c->next[2*j+b];
    }
    *id=rank;return 0;
}
int uc_decode(const UC *c,uint64_t id,uint64_t *path) {
    unsigned j=3*c->m/8,r=c->R,n=0;uint64_t word=1;
    if(id>=c->counts[(size_t)r*c->J+j])return -1;
    while(id) {
        id--;uint64_t left=count0(c,j,r);unsigned b=id>=left;
        if(b)id-=left;
        unsigned k=c->cost[2*j+b];if(k>r||++n>c->p)return -2;
        r-=k;j=c->next[2*j+b];word=(word<<1)|b;
    }
    *path=word;return 0;
}
/* Precondition: state is ID 0 initially or an ID produced by this kernel.
   Quantized states and budgets are immutable model data, not hidden state. */
int uc_advance(const UC *c,uint8_t *state,unsigned nb,unsigned b,uint64_t out[4]) {
    if(c->p<4||c->p>30||!nb||nb>4||b>1)return -1;
    uint64_t path;if(uc_decode(c,load(state,nb),&path))return -2;
    Cursor old=root();old.path=path;
    for(uint64_t v=path;v>1;v>>=1)old.depth++;
    if(old.depth>20)return -2;
    Cursor z=update(old,b,c->p);uint64_t id;
    if(uc_encode(c,z.path,&id))return -3;
    if(id>=((uint64_t)1<<(8*nb)))return -4;
    save(state,nb,id);output(z,c->p,out);return 0;
}
int uc_output(const UC *c,const uint8_t *state,unsigned nb,uint64_t out[4]) {
    if(c->p<4||c->p>30||!nb||nb>4)return -1;
    uint64_t path;if(uc_decode(c,load(state,nb),&path))return -2;
    Cursor z;if(from_path(c->p,path,&z))return -2;
    output(z,c->p,out);return 0;
}
