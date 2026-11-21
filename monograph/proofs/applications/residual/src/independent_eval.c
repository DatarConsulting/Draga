#include <math.h>
#include <stddef.h>
void eval(const float *u,const float *v,int rank,float *out){
 for(int i=0;i<128*128;i++)out[i]=0;
 for(int a=0;a<rank;a++)for(int i=0;i<128;i++)for(int j=0;j<128;j++)
  out[i*128+j]=fmaf(u[a*128+i],v[a*128+j],out[i*128+j]);
}
