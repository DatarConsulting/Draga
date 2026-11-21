#include <stdint.h>
#include <stddef.h>
typedef unsigned __int128 u128;
typedef __int128 i128;
typedef struct { uint64_t a,b,c,d; unsigned depth; uint64_t path; } Cursor;
/* Fixed source only. p and bytes are immutable caller parameters. No table,
   allocation, mutable global, or external clock. p in [4,30], depth <=20. */
static Cursor root(void) { Cursor z={1,0,0,1,0,1}; return z; }
static Cursor append(Cursor z,unsigned b) {
    uint64_t a=z.a,c=z.c;
    if(b==0) { z.a=2*a+z.b; z.b=3*a+6*z.b; z.c=2*c+z.d; z.d=3*c+6*z.d; }
    else { z.a=6*a+3*z.b; z.b=a+2*z.b; z.c=6*c+3*z.d; z.d=c+2*z.d; }
    z.path=(z.path<<1)|b; z.depth++; return z;
}
static Cursor pop(Cursor z) {
    uint64_t a=z.a,b=z.b,c=z.c,d=z.d;
    if((z.path&1)==0) {
        z.a=(uint64_t)(((i128)6*a-b)/9); z.b=(uint64_t)(((i128)2*b-3*a)/9);
        z.c=(uint64_t)(((i128)6*c-d)/9); z.d=(uint64_t)(((i128)2*d-3*c)/9);
    } else {
        z.a=(uint64_t)(((i128)2*a-3*b)/9); z.b=(uint64_t)(((i128)6*b-a)/9);
        z.c=(uint64_t)(((i128)2*c-3*d)/9); z.d=(uint64_t)(((i128)6*d-c)/9);
    }
    z.path>>=1;z.depth--;return z;
}
static int leaf(Cursor z,unsigned p) {
    u128 det=(u128)z.a*z.d-(u128)z.b*z.c;
    u128 dl=3*((u128)z.a+z.c)+5*((u128)z.b+z.d);
    u128 dr=5*((u128)z.a+z.c)+3*((u128)z.b+z.d);
    return (det*528*((uint64_t)1<<p))<=49*dl*dr;
}
static int next(Cursor *z,unsigned p) {
    if(!leaf(*z,p)) { *z=append(*z,0);return 1; }
    while(z->depth) {
        unsigned b=(unsigned)(z->path&1);*z=pop(*z);
        if(!b) { *z=append(*z,1);return 1; }
    }
    return 0;
}
static int decode(unsigned p,uint64_t index,Cursor *result) {
    Cursor z=root();uint64_t i=0;
    do { if(i++==index) { *result=z;return 0; } } while(next(&z,p));
    return -1;
}
static int rank_path(unsigned p,uint64_t path,uint64_t *result) {
    Cursor z=root();uint64_t i=0;
    do { if(z.path==path) { *result=i;return 0; } i++; } while(next(&z,p));
    return -1;
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
int cf_count(unsigned p,uint64_t *n,unsigned *depth) {
    if(p<4||p>30)return -1;Cursor z=root();*n=0;*depth=0;
    do {(*n)++;if(z.depth>*depth)*depth=z.depth;}while(next(&z,p));return 0;
}
int cf_rank_step(unsigned p,uint8_t *state,unsigned nb,unsigned b) {
    if(p<4||p>30||nb==0||nb>4||b>1)return -1;
    Cursor z; if(decode(p,load(state,nb),&z))return -2;
    z=update(z,b,p);uint64_t i;if(rank_path(p,z.path,&i))return -3;
    if(i>=((uint64_t)1<<(8*nb)))return -4;save(state,nb,i);return 0;
}
static int from_path(unsigned p,uint64_t path,Cursor *result) {
    if(path==0)return -1;unsigned depth=0;for(uint64_t t=path;t>1;t>>=1)depth++;
    if(depth>20)return -1;Cursor z=root();
    for(unsigned i=depth;i>0;i--) {if(leaf(z,p))return -1;z=append(z,(path>>(i-1))&1);}
    *result=z;return 0;
}
int cf_path_step(unsigned p,uint8_t *state,unsigned nb,unsigned b) {
    if(p<4||p>30||nb==0||nb>4||b>1)return -1;
    Cursor z;if(from_path(p,load(state,nb),&z))return -2;z=update(z,b,p);
    if(z.path>=((uint64_t)1<<(8*nb)))return -4;save(state,nb,z.path);return 0;
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
int cf_rank_output(unsigned p,const uint8_t *state,unsigned nb,uint64_t out[4]) {
    if(p<4||p>30||nb==0||nb>4)return -1;Cursor z;
    if(decode(p,load(state,nb),&z))return -2;output(z,p,out);return 0;
}
int cf_path_output(unsigned p,const uint8_t *state,unsigned nb,uint64_t out[4]) {
    if(p<4||p>30||nb==0||nb>4)return -1;Cursor z;
    if(from_path(p,load(state,nb),&z))return -2;output(z,p,out);return 0;
}
