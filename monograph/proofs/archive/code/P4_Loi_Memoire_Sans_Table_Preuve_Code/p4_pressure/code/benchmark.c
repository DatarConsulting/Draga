#define _POSIX_C_SOURCE 200809L
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include "compact_filter.c"
#include "../reference/packed_predictor.c"
static double now(void) {struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+1e-9*t.tv_nsec;}
static uint64_t rng(uint64_t *z) {*z^=*z<<13;*z^=*z>>7;*z^=*z<<17;return *z;}
int main(int argc,char **argv) {
 if(argc!=2)return 1;FILE *f=fopen(argv[1],"rb");if(!f)return 2;
 fseek(f,0,SEEK_END);long n=ftell(f);rewind(f);uint8_t *m=malloc((size_t)n);
 if(!m||fread(m,1,(size_t)n,f)!=(size_t)n)return 3;fclose(f);
 if(!pc_valid(m,(size_t)n-2))return 4;
 const unsigned p=19,nb=2,steps=20000;
 for(unsigned repeat=0;repeat<3;repeat++)for(unsigned mode=0;mode<3;mode++) {
  uint8_t state[2]={0,0};if(mode==2)state[0]=1;
  uint64_t seed=0x123456789abcdefull,out[4],checksum=0;double start=now();
  for(unsigned i=0;i<steps;i++) {
   unsigned b=rng(&seed)&1;int rc;
   if(mode==0){rc=pc_step(m,state,b);rc|=pc_output_stream(m,state,out);}
   else if(mode==1){rc=cf_rank_step(p,state,nb,b);rc|=cf_rank_output(p,state,nb,out);}
   else {rc=cf_path_step(p,state,nb,b);rc|=cf_path_output(p,state,nb,out);}
   if(rc)return 5;checksum^=out[0]+state[0]+(i*0x9e3779b9u);
  }
  double elapsed=now()-start;
  printf("%u,%u,%u,%.9f,%.3f,%llu\n",repeat,mode,steps,elapsed,elapsed*1e9/steps,(unsigned long long)checksum);
 }
 free(m);return 0;
}
