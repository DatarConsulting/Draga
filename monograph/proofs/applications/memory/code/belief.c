#include <stdint.h>
/* Exact normalized posterior on a grid. Requires 0<=j<=G<=2^22 and row sum 2^24. */
uint32_t belief_step(uint32_t j,int a,uint32_t G,int K,const uint32_t *A) {
    int stride=2*K+1;
    __uint128_t n=(__uint128_t)j*A[2*a]+(__uint128_t)(G-j)*A[stride+2*a];
    __uint128_t d=n+(__uint128_t)j*A[2*a+1]+(__uint128_t)(G-j)*A[stride+2*a+1];
    __uint128_t num=n*G,q=num/d,r=num%d;
    if(2*r>d || (2*r==d && (q&1))) q++;
    return (uint32_t)q;
}
int belief_step_packed(uint8_t *state,int nbytes,int a,uint32_t G,int K,const uint32_t *A) {
    if(nbytes<1 || nbytes>4)return -1;
    uint32_t j=0;
    for(int i=0;i<nbytes;i++)j|=((uint32_t)state[i])<<(8*i);
    if(j>G || a<0 || a>=K)return -2;
    j=belief_step(j,a,G,K,A);
    if(nbytes<4 && j>=((uint32_t)1<<(8*nbytes)))return -3;
    for(int i=0;i<nbytes;i++)state[i]=(j>>(8*i))&255;
    return 0;
}
uint64_t belief_probabilities(const uint8_t *state,int nbytes,uint32_t G,int K,
                              const uint32_t *A,uint64_t *numerators) {
    uint32_t j=0;int stride=2*K+1;
    for(int i=0;i<nbytes;i++)j|=((uint32_t)state[i])<<(8*i);
    for(int a=0;a<K;a++)numerators[a]=(uint64_t)j*(A[2*a]+A[2*a+1])+(uint64_t)(G-j)*(A[stride+2*a]+A[stride+2*a+1]);
    numerators[K]=(uint64_t)j*A[2*K]+(uint64_t)(G-j)*A[stride+2*K];
    return (uint64_t)G*((uint64_t)1<<24);
}
/* Local smoothing with lambda=2^-s; verified callers bound all uint64 sums. */
uint64_t belief_probabilities_smoothed(const uint8_t *state,int nbytes,uint32_t G,
                                      int K,const uint32_t *A,int s,uint64_t *nums) {
    uint64_t den=belief_probabilities(state,nbytes,G,K,A,nums),L=(uint64_t)1<<s;
    for(int a=0;a<=K;a++)nums[a]=(L-1)*(2*K)*nums[a]+den*(a==K?K:1);
    return den*L*(2*K);
}
