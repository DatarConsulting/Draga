#include <stdint.h>
#include <stddef.h>
/* All mutable inter-call memory is the caller's packed byte array.
   a,m are immutable model parameters. No external clock or history is read. */
static int width(int a,int m,int t) {
    int v=(a+1)*m+1-a*t;
    return v<t?v:t;
}
static int config(int a,int m,int *L,int *T,uint64_t *N) {
    if(a<1 || a>30 || m<a+2 || m>30 || (a+1)*m+2>62)return -1;
    *L=(a+1)*m+2;*T=(*L-2+a)/a;*N=1;
    for(int t=0;t<*T;t++)*N+=(uint64_t)1<<width(a,m,t);
    return 0;
}
static uint64_t offset(int a,int m,int t) {
    uint64_t v=0;for(int j=0;j<t;j++)v+=(uint64_t)1<<width(a,m,j);return v;
}
static uint64_t read_state(const uint8_t *p,int n) {
    uint64_t v=0;for(int i=0;i<n;i++)v|=(uint64_t)p[i]<<(8*i);return v;
}
static int decode(int a,int m,int T,uint64_t N,uint64_t id,int *t,uint64_t *b) {
    if(id>=N)return -1;
    if(id==N-1){*t=T;*b=0;return 0;}
    for(*t=0;*t<T;(*t)++){
        uint64_t sz=(uint64_t)1<<width(a,m,*t);
        if(id<sz){*b=id;return 0;}id-=sz;
    }
    return -1;
}
int p4_step(uint8_t *p,int nbytes,int a,int m,int token) {
    int L,T,t;uint64_t N,b;
    if(config(a,m,&L,&T,&N) || nbytes!=(m+9)/8 || token<0 || token>1)return -1;
    if(decode(a,m,T,N,read_state(p,nbytes),&t,&b))return -2;
    uint64_t next=N-1;
    if(t<T-1){
        int oldn=width(a,m,t),newn=width(a,m,t+1);
        uint64_t nb=((uint64_t)token<<(newn-1))+(b>>(oldn-newn+1));
        next=offset(a,m,t+1)+nb;
    }
    for(int i=0;i<nbytes;i++)p[i]=(uint8_t)(next>>(8*i));
    return 0;
}
int p4_output(const uint8_t *p,int nbytes,int a,int m,uint64_t *numerator,int *power) {
    int L,T,t;uint64_t N,b;
    if(config(a,m,&L,&T,&N) || nbytes!=(m+9)/8)return -1;
    if(decode(a,m,T,N,read_state(p,nbytes),&t,&b))return -2;
    *power=L;
    if(t==T)*numerator=1;
    else if(t<=m)*numerator=b<<(L-(a+1)*t);
    else *numerator=2*b+1;
    return 0;
}
#ifdef P4_SANITIZER_MAIN
int main(void){
    for(int a=1;a<=6;a++)for(int m=a+2;m<=12 && (a+1)*m+2<=62;m++){
        int L,T;uint64_t N;config(a,m,&L,&T,&N);int nb=(m+9)/8;
        for(uint64_t s=0;s<N;s++)for(int token=0;token<2;token++){
            uint8_t b[8]={0};for(int j=0;j<nb;j++)b[j]=(uint8_t)(s>>(8*j));
            uint64_t out;int p;
            if(p4_output(b,nb,a,m,&out,&p) || p4_step(b,nb,a,m,token))return 1;
            if(p4_output(b,nb,a,m,&out,&p))return 2;
        }
    }
    return 0;
}
#endif
