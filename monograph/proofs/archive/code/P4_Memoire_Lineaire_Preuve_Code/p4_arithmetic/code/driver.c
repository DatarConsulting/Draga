#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "arithmetic_filter.c"
static double now(void){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+1e-9*t.tv_nsec;}
static uint64_t rng(uint64_t*s){*s^=*s<<13;*s^=*s>>7;*s^=*s<<17;return *s;}
int main(int argc,char**argv){
    if(argc!=3)return 1;
    if(argv[1][0]=='v'){
        FILE*f=fopen(argv[2],"r");if(!f)return 2;
        unsigned p,B,nb,n;if(fscanf(f,"%u %u %u %u",&p,&B,&nb,&n)!=4)return 3;
        for(unsigned k=0;k<n;k++){
            unsigned long long id,dest[2];if(fscanf(f,"%llu %llu %llu",&id,&dest[0],&dest[1])!=3)return 4;
            uint64_t path,back,out[4];uint8_t state[4];
            if(ac_decode(p,id,&path)||ac_encode(p,path,&back)||back!=id)return 5;
            save(state,nb,id);if(ac_output(p,state,nb,out))return 6;
            for(unsigned b=0;b<2;b++){save(state,nb,id);if(ac_advance(p,state,nb,b,out)||load(state,nb)!=dest[b])return 7;}
        }
        fclose(f);printf("p=%u states=%u transitions=%u ASan/UBSan pass\n",p,n,2*n);return 0;
    }
    unsigned p=(unsigned)strtoul(argv[2],NULL,10),nb=(ac_bits(p)+7)/8,steps=100000;
    for(unsigned rep=0;rep<5;rep++){
        uint8_t state[4]={0};uint64_t out[4],seed=0x123456789abcdefull,checksum=0;double start=now();
        for(unsigned i=0;i<steps;i++){
            if(ac_advance(p,state,nb,rng(&seed)&1,out))return 8;checksum^=out[0]+i;
        }
        printf("%u,%u,%u,%.3f,%llu\n",p,rep,steps,(now()-start)*1e9/steps,(unsigned long long)checksum);
    }
    return 0;
}
