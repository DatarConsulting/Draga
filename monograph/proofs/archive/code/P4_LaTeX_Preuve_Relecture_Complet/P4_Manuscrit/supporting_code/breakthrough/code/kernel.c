#include <stdint.h>
#include <stddef.h>
#include <limits.h>

/* Packed persistent state only. k<=4, total bits<=64, 3<=bits[i]<=32.
   Intermediate signed 128-bit accumulation; caller verifies its exact bound.
   Rounding: nearest, ties to even. No floating point in the recurrence. */
static int64_t nearest_shift(__int128 v,int sh) {
    if(sh<=0) return (int64_t)(v*(((__int128)1)<<(-sh)));
    int sign=v<0?-1:1;
    unsigned __int128 a=v<0?(unsigned __int128)(-v):(unsigned __int128)v;
    unsigned __int128 q=a>>sh;
    unsigned __int128 r=a- (q<<sh);
    unsigned __int128 half=((unsigned __int128)1)<<(sh-1);
    if(r>half || (r==half && (q&1))) q++;
    return sign*(int64_t)q;
}
int step(uint8_t *state,int nbytes,int k,const int32_t *a,
         const int32_t *bits,const int32_t *shift,int token,int qbits) {
    if(k<1||k>4||nbytes<1||nbytes>8||token<0||token>1)return -1;
    uint64_t packed=0;
    for(int z=0;z<nbytes;z++)packed|=((uint64_t)state[z])<<(8*z);
    int64_t old[4],next[4];int off=0,minshift=shift[0]-qbits;
    for(int i=0;i<k;i++){
      if(bits[i]<3||bits[i]>32)return -2;
      uint64_t mask=(((uint64_t)1)<<bits[i])-1;
      uint64_t raw=(packed>>off)&mask;
      old[i]=(raw&(((uint64_t)1)<<(bits[i]-1)))?(int64_t)raw-(((int64_t)1)<<bits[i]):(int64_t)raw;
      off+=bits[i];if(shift[i]-qbits<minshift)minshift=shift[i]-qbits;
    }
    if(off>8*nbytes)return -3;
    for(int j=0;j<k;j++){
      __int128 sum=0;
      for(int i=0;i<k;i++){
        int d=shift[i]-qbits-minshift;if(d<0||d>60)return -4;
        __int128 v=(__int128)old[i]*a[token*k*k+i*k+j];
        sum+=v*(((__int128)1)<<d);
      }
      int sh=shift[j]-minshift;if(sh<0||sh>100)return -5;
      next[j]=nearest_shift(sum,sh);
      if(next[j]<-(((int64_t)1)<<(bits[j]-1))||next[j]>=(((int64_t)1)<<(bits[j]-1)))return 1;
    }
    uint64_t out=0;off=0;
    for(int i=0;i<k;i++){
      uint64_t mask=(((uint64_t)1)<<bits[i])-1;
      out|=(((uint64_t)next[i])&mask)<<off;off+=bits[i];
    }
    for(int z=0;z<nbytes;z++)state[z]=(uint8_t)(out>>(8*z));
    return 0;
}
/* Exact dyadic output numerator * 2^exponent; both are output temporaries. */
int output_value(const uint8_t *state,int nbytes,int k,const int32_t *beta,
                 const int32_t *bits,const int32_t *shift,int qbits,
                 int64_t *numerator,int *exponent) {
    uint64_t packed=0;int off=0,base=shift[0];__int128 sum=0;
    for(int i=0;i<k;i++)if(shift[i]<base)base=shift[i];
    for(int z=0;z<nbytes;z++)packed|=((uint64_t)state[z])<<(8*z);
    for(int i=0;i<k;i++){
        uint64_t mask=(((uint64_t)1)<<bits[i])-1,raw=(packed>>off)&mask;
        int64_t value=(raw&(((uint64_t)1)<<(bits[i]-1)))?(int64_t)raw-(((int64_t)1)<<bits[i]):(int64_t)raw;
        sum+=(__int128)value*beta[i]*(((__int128)1)<<(shift[i]-base));off+=bits[i];
    }
    if(sum>INT64_MAX || sum<INT64_MIN)return 1;
    *numerator=(int64_t)sum;*exponent=base-qbits;return 0;
}
/* One-byte identifier for a fully enumerated machine with <=256 states. */
int table_step(uint8_t *state,int token,int count,const uint8_t *next) {
    if(token<0||token>1||*state>=count)return -1;
    *state=next[2*(*state)+token];return 0;
}
