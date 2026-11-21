#include <stdint.h>
#include <stddef.h>
typedef unsigned __int128 u128;
typedef struct { uint64_t a,b,c,d; unsigned depth; uint64_t path; } Cursor;
/* Fixed source only. p and bytes are immutable caller parameters. No table,
   allocation, mutable global, or external clock. p in [4,22], depth <=12; width proof inherited for p<=30. */
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
/* Injective relabelling, not the old preorder rank. For p=19,21, exactly
   four deepest paths are moved into four absent heap codes. No tree data. */
static uint64_t unfold(unsigned p,uint64_t id) {
    if(p!=19 && p!=21) return id;
    unsigned B=(p+1)/2;
    uint64_t H=(uint64_t)1<<(B-1), v=(H-1)/3;
    if(id<H || id>=H+4) return id;
    uint64_t offset=id-H;
    uint64_t parent=H+((offset>>1)?H-1-v:v);
    return (parent<<1)|(offset&1);
}
static uint64_t fold(unsigned p,uint64_t path) {
    if(p!=19 && p!=21) return path;
    unsigned B=(p+1)/2;
    uint64_t H=(uint64_t)1<<(B-1), v=(H-1)/3;
    if((path>>1)==H+v) return H+(path&1);
    if((path>>1)==2*H-1-v) return H+2+(path&1);
    return path;
}
/* Canonical, reachable IDs only. Bytes above the advertised B-bit budget
   must be zero; from_path checks actual validity. Unused IDs are rejected. */
static int decode_id(unsigned p,uint64_t id,Cursor *z) {
    if(from_path(p,unfold(p,id),z)) return -2;
    if(fold(p,z->path)!=id) return -2;
    return 0;
}
int df_step(unsigned p,uint8_t *state,unsigned nb,unsigned b) {
    if(p<4||p>22||nb==0||nb>4||b>1) return -1;
    Cursor z;
    if(decode_id(p,load(state,nb),&z)) return -2;
    z=update(z,b,p);
    uint64_t id=fold(p,z.path);
    if(id>=((uint64_t)1<<(8*nb))) return -4;
    save(state,nb,id);
    return 0;
}
int df_output(unsigned p,const uint8_t *state,unsigned nb,uint64_t out[4]) {
    if(p<4||p>22||nb==0||nb>4) return -1;
    Cursor z;
    if(decode_id(p,load(state,nb),&z)) return -2;
    output(z,p,out);
    return 0;
}
/* Same transition and exact readout, sharing the reconstructed matrix. */
int df_step_output(unsigned p,uint8_t *state,unsigned nb,unsigned b,uint64_t out[4]) {
    if(p<4||p>22||nb==0||nb>4||b>1) return -1;
    Cursor z;
    if(decode_id(p,load(state,nb),&z)) return -2;
    z=update(z,b,p);
    uint64_t id=fold(p,z.path);
    if(id>=((uint64_t)1<<(8*nb))) return -4;
    save(state,nb,id);
    output(z,p,out);
    return 0;
}
/* Fast production API. Precondition: state is initial ID 1 or a reachable
   ID produced by this kernel, with nb sufficient. The caller owns p,nb.
   It does not revalidate the old suffix; it directly constructs the new one.
   No old matrix is reconstructed. All reachable transitions are certified. */
int df_advance(unsigned p,uint8_t *state,unsigned nb,unsigned b,uint64_t out[4]) {
    if(p<4||p>22||nb==0||nb>4||b>1) return -1;
    uint64_t path=unfold(p,load(state,nb));
    if(!path) return -2;
    unsigned depth=0;for(uint64_t t=path;t>1;t>>=1)depth++;
    if(depth>20) return -2;
    Cursor old=root();old.path=path;old.depth=depth;
    Cursor z=update(old,b,p);
    uint64_t id=fold(p,z.path);
    if(id>=((uint64_t)1<<(8*nb))) return -4;
    save(state,nb,id);output(z,p,out);
    return 0;
}
