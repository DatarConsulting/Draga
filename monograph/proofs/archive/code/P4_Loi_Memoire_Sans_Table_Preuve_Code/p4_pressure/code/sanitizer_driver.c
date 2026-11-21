#include <stdio.h>
#include "compact_filter.c"
int main(void) {
 for(unsigned p=4;p<=22;p++) {
  uint64_t n;unsigned depth;if(cf_count(p,&n,&depth))return 1;
  Cursor z=root();uint64_t i=0;unsigned nb=1;while(n>((uint64_t)1<<(8*nb)))nb++;
  unsigned nbp=(depth+8)/8;
  do {
   for(unsigned mode=0;mode<2;mode++)for(unsigned b=0;b<2;b++) {
    uint8_t state[4]={0};uint64_t out[4];unsigned bytes=mode?nbp:nb;
    save(state,bytes,mode?z.path:i);
    int rc=mode?cf_path_output(p,state,bytes,out):cf_rank_output(p,state,bytes,out);
    rc|=mode?cf_path_step(p,state,bytes,b):cf_rank_step(p,state,bytes,b);
    if(rc)return 2;
   }
   i++;
  }while(next(&z,p));
  if(i!=n)return 3;
 }
 puts("Address and undefined behavior checks passed for all 19 finite graphs.");
 return 0;
}
