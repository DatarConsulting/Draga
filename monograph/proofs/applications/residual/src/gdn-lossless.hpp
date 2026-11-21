#pragma once
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <stdexcept>

// Exact IEEE-754 bit representation: sign/mantissa (24 bits), and rowwise
// bit-packed exponent offsets. A raw escape handles incompressible heads.
// No floating-point arithmetic, rounding, future query, or calibration.
namespace gdn_lossless {
constexpr size_t elements = 128*128;
constexpr size_t raw_bytes = elements*4;
constexpr size_t max_encoded = raw_bytes+1;
constexpr size_t scratch_bytes = raw_bytes+2*128+1;
inline size_t encode(const float * x, uint8_t * out) {
    bool zero=true;
    for(size_t i=0;i<elements;++i) {uint32_t u;std::memcpy(&u,x+i,4);zero &= u==0;}
    if(zero) return 0; // implicit positive zero, preserves negative-zero bits
    out[0]=1;size_t pos=1;
    for(size_t row=0;row<128;++row) {
        uint32_t words[128];unsigned lo=255,hi=0;
        for(size_t i=0;i<128;++i) {
            std::memcpy(words+i,x+row*128+i,4);
            const unsigned e=(words[i]>>23)&255;lo=std::min(lo,e);hi=std::max(hi,e);
        }
        unsigned bits=0;while((1u<<bits)<=hi-lo) ++bits;
        out[pos++]=uint8_t(lo);out[pos++]=uint8_t(bits);
        for(uint32_t u:words) {
            const uint32_t m=(u&0x7fffff)|((u>>8)&0x800000);
            out[pos++]=uint8_t(m);out[pos++]=uint8_t(m>>8);out[pos++]=uint8_t(m>>16);
        }
        uint32_t acc=0;unsigned have=0;
        for(uint32_t u:words) {
            acc|=(((u>>23)&255)-lo)<<have;have+=bits;
            while(have>=8) {out[pos++]=uint8_t(acc);acc>>=8;have-=8;}
        }
        if(have) throw std::logic_error("unaligned exponent row");
    }
    if(pos>=max_encoded) {out[0]=0;std::memcpy(out+1,x,raw_bytes);return max_encoded;}
    return pos;
}
inline void decode(const uint8_t * in,size_t n,float * x) {
    if(n==0) {std::memset(x,0,raw_bytes);return;}
    if(n>max_encoded) throw std::runtime_error("oversized lossless head");
    if(in[0]==0) {
        if(n!=max_encoded) throw std::runtime_error("invalid raw head size");
        std::memcpy(x,in+1,raw_bytes);return;
    }
    if(in[0]!=1) throw std::runtime_error("invalid lossless tag");
    size_t pos=1;
    for(size_t row=0;row<128;++row) {
        if(n-pos<2) throw std::runtime_error("truncated lossless row");
        const unsigned lo=in[pos++],bits=in[pos++];
        if(bits>8 || n-pos<384+16*bits) throw std::runtime_error("invalid lossless row size");
        const uint8_t * mant=in+pos;pos+=384;
        uint32_t acc=0;unsigned have=0;const uint32_t mask=(1u<<bits)-1;
        for(size_t i=0;i<128;++i) {
            if(have<bits) {acc|=uint32_t(in[pos++])<<have;have+=8;}
            const unsigned exp=lo+(acc&mask);acc>>=bits;have-=bits;
            if(exp>255) throw std::runtime_error("invalid lossless exponent");
            const uint32_t m=uint32_t(mant[3*i])|(uint32_t(mant[3*i+1])<<8)|(uint32_t(mant[3*i+2])<<16);
            const uint32_t u=(m&0x7fffff)|((m&0x800000)<<8)|(exp<<23);
            std::memcpy(x+row*128+i,&u,4);
        }
    }
    if(pos!=n) throw std::runtime_error("trailing lossless bytes");
}
}
