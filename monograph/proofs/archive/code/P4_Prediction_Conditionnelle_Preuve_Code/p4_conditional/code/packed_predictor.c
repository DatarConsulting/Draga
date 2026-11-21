#include <stdint.h>
#include <stddef.h>
/* Immutable model: PCV1 header (16 bytes), nodes (j:u32,next0:u16,next1:u16).
   The ONLY persistent mutable data is state[nb]. No allocation or hidden clock. */
static uint16_t u16(const uint8_t *p) { return (uint16_t)p[0] | ((uint16_t)p[1]<<8); }
static uint32_t u32(const uint8_t *p) { return (uint32_t)p[0] | ((uint32_t)p[1]<<8) | ((uint32_t)p[2]<<16) | ((uint32_t)p[3]<<24); }
static uint32_t sid(const uint8_t *s, unsigned nb) { uint32_t z=0; for(unsigned i=0;i<nb;i++) z|=(uint32_t)s[i]<<(8*i); return z; }
int pc_valid(const uint8_t *m, size_t size) {
    if(size<16 || m[0]!='P'||m[1]!='C'||m[2]!='V'||m[3]!='1') return 0;
    unsigned p=m[4],nb=m[5]; uint32_t n=u16(m+6),bits=u32(m+8);
    if(p<4||p>22||n==0||bits==0||bits>16||nb!=(bits+7)/8) return 0;
    if(size!=16+8*(size_t)n || n>(1u<<bits) || u32(m+12)>=n) return 0;
    for(uint32_t i=0;i<n;i++) {
        const uint8_t *v=m+16+8*i;
        if(u32(v)>(1u<<p) || u16(v+4)>=n || u16(v+6)>=n) return 0;
    }
    return 1;
}
int pc_init(const uint8_t *m,uint8_t *state) {
    uint32_t z=u32(m+12); for(unsigned i=0;i<m[5];i++) state[i]=(uint8_t)(z>>(8*i)); return 0;
}
int pc_step(const uint8_t *m,uint8_t *state,unsigned b) {
    uint32_t z=sid(state,m[5]); if(b>1||z>=u16(m+6)) return -1;
    z=u16(m+16+8*z+4+2*b);
    for(unsigned i=0;i<m[5];i++) state[i]=(uint8_t)(z>>(8*i));
    return 0;
}
int pc_output(const uint8_t *m,const uint8_t *state,uint64_t *out) {
    uint32_t z=sid(state,m[5]); if(z>=u16(m+6)) return -1;
    uint64_t D=1ull<<m[4],j=u32(m+16+8*z);
    out[0]=3*D-2*j; out[1]=D+2*j; out[2]=4*D; out[3]=8*D;
    return 0;
}
/* Stationary binary stream, with no termination: exact dyadic probabilities. */
int pc_output_stream(const uint8_t *m,const uint8_t *state,uint64_t *out) {
    uint32_t z=sid(state,m[5]); if(z>=u16(m+6)) return -1;
    uint64_t D=1ull<<m[4],j=u32(m+16+8*z);
    out[0]=3*D-2*j; out[1]=D+2*j; out[2]=4*D;
    return 0;
}
