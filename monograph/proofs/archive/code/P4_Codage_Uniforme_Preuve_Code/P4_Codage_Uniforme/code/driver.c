#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "uniform_filter.c"
int df_advance(unsigned,uint8_t*,unsigned,unsigned,uint64_t*);
static double now(void){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+1e-9*t.tv_nsec;}
static uint64_t rng(uint64_t *s){*s^=*s<<13;*s^=*s>>7;*s^=*s<<17;return *s;}
int main(int argc,char **argv){
    if(argc!=3)return 1;FILE*f=fopen(argv[2],"rb");if(!f)return 2;
    uint32_t h[5];if(fread(h,4,5,f)!=5)return 3;
    UC c={h[0],h[1],h[2],h[3],NULL,NULL,NULL};
    uint32_t *next=malloc(8*c.J),*cost=malloc(8*c.J);
    uint64_t *counts=malloc(8*(size_t)c.J*(c.R+1));
    if(!next||!cost||!counts)return 4;
    if(fread(next,4,2*c.J,f)!=2*c.J||fread(cost,4,2*c.J,f)!=2*c.J)return 5;
    fclose(f);c.next=next;c.cost=cost;c.counts=counts;
    if(uc_build(c.p,c.J,c.R,next,cost,counts))return 6;
    uint64_t total=counts[(size_t)c.R*c.J+3*c.m/8];
    unsigned nb=1;while(total>((uint64_t)1<<(8*nb)))nb++;
    if(argv[1][0]=='v'){
        uint64_t accepted=0,edges=0;
        for(uint64_t i=0;i<total;i++){
            uint64_t path,id;if(uc_decode(&c,i,&path)||uc_encode(&c,path,&id)||id!=i)return 7;
            uint8_t s[4]={0};uint64_t out[4];save(s,nb,i);
            int rc=uc_output(&c,s,nb,out);if(rc==-2)continue;if(rc)return 8;
            accepted++;
            for(unsigned b=0;b<2;b++){
                save(s,nb,i);if(uc_advance(&c,s,nb,b,out))return 9;edges++;
            }
        }
        if(accepted!=h[4])return 10;
        printf("p=%u: %llu reachable, %llu transitions, %llu language IDs; ASan/UBSan passed.\n",c.p,(unsigned long long)accepted,(unsigned long long)edges,(unsigned long long)total);
    }else{
        unsigned modes=c.p==23?2:1,steps=100000;
        for(unsigned rep=0;rep<5;rep++)for(unsigned k=0;k<modes;k++){
            unsigned mode=(k+rep)%modes;uint8_t s[4]={0};if(mode)s[0]=1;
            uint64_t seed=0x123456789abcdefull,out[4],checksum=0;double start=now();
            for(unsigned i=0;i<steps;i++){
                unsigned b=rng(&seed)&1;
                int rc=mode?df_advance(23,s,2,b,out):uc_advance(&c,s,nb,b,out);
                if(rc)return 11;checksum^=out[0]+i;
            }
            double elapsed=now()-start;
            printf("%u,%u,%u,%u,%.3f,%llu\n",c.p,rep,mode,steps,elapsed*1e9/steps,(unsigned long long)checksum);
        }
    }
    free(next);free(cost);free(counts);return 0;
}
