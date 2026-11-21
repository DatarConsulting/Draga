#include <stdint.h>
#include <stdio.h>
#include <string.h>
int cf_count(unsigned,uint64_t*,unsigned*);
int df_step(unsigned,uint8_t*,unsigned,unsigned);
int df_output(unsigned,const uint8_t*,unsigned,uint64_t*);
int df_step_output(unsigned,uint8_t*,unsigned,unsigned,uint64_t*);
int df_advance(unsigned,uint8_t*,unsigned,unsigned,uint64_t*);
int main(void) {
    uint64_t checked=0;
    for(unsigned p=4;p<=22;p++) {
        uint64_t n;unsigned depth;
        if(cf_count(p,&n,&depth))return 1;
        unsigned B=0;while(((uint64_t)1<<B)<n)B++;
        unsigned nb=(B+7)/8;uint64_t accepted=0;
        for(uint64_t id=0;id<((uint64_t)1<<B);id++) {
            uint8_t state[4]={0};uint64_t out[4];
            for(unsigned k=0;k<nb;k++)state[k]=(uint8_t)(id>>(8*k));
            int rc=df_output(p,state,nb,out);
            if(rc==-2)continue;
            if(rc)return 2;
            accepted++;
            for(unsigned b=0;b<2;b++) {
                uint8_t a[4],c[4];uint64_t oa[4],oc[4];
                memcpy(a,state,4);memcpy(c,state,4);
                if(df_step(p,a,nb,b)||df_output(p,a,nb,oa)||df_step_output(p,c,nb,b,oc))return 3;
                if(memcmp(a,c,4)||memcmp(oa,oc,sizeof(oa)))return 4;
                memcpy(c,state,4);
                if(df_advance(p,c,nb,b,oc)||memcmp(a,c,4)||memcmp(oa,oc,sizeof(oa)))return 6;
                checked+=3;
            }
        }
        if(accepted!=n)return 5;
    }
    printf("Address/UB checks passed: %llu transitions.\n",(unsigned long long)checked);
    return 0;
}
