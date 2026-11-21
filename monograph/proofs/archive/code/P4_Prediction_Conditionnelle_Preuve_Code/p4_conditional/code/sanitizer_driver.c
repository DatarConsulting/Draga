#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "packed_predictor.c"
int main(int argc,char **argv) {
    if(argc<2) return 2;
    for(int a=1;a<argc;a++) {
        FILE *f=fopen(argv[a],"rb"); if(!f) return 3;
        if(fseek(f,0,SEEK_END)) return 4; long n=ftell(f); if(n<16) return 5; rewind(f);
        uint8_t *m=malloc((size_t)n); if(!m) return 6;
        if(fread(m,1,(size_t)n,f)!=(size_t)n) return 7; fclose(f);
        unsigned nb=m[5]; if(!pc_valid(m,(size_t)n-nb)) return 8;
        uint8_t state[2]; uint64_t out[4];
        for(uint32_t i=0;i<u16(m+6);i++) {
            for(unsigned b=0;b<2;b++) {
                for(unsigned k=0;k<nb;k++) state[k]=(uint8_t)(i>>(8*k));
                if(pc_output(m,state,out)||pc_output_stream(m,state,out)||pc_step(m,state,b)) return 9;
            }
        }
        free(m);
    }
    puts("All packed states and transitions passed address/undefined checks.");
    return 0;
}
