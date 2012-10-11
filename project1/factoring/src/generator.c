#include <stdio.h>
#include <stdlib.h>

#include <gmp.h>

int main(int argc,char* argv[])
{
    
    if(argc!=3) {
        printf("Wrong number of arguments: %d",argc-1);
        return 1;
    }
    int bits;
    int N;
    bits = atoi(argv[1]);
    N = atoi(argv[2]);
    
    mpz_t number;
    gmp_randstate_t rnd;
    gmp_randinit_default(rnd);

    mpz_init(number);

    int i;
    for(i=0;i<N;++i)
    {
        mpz_urandomb(number,rnd,bits);
        gmp_printf("%Zd\n",number);
    }

    return 0;
    
}
