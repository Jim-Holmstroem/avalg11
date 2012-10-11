#include <stdio.h>
#include <stdlib.h>

#include <gmp.h>

int main()
{

    FILE *in = NULL;
    FILE *out = NULL;

    char name_in[128];
    char name_out[128];

    gmp_randstate_t rnd;
    gmp_randinit_default(rnd);
    mpz_t a;
    mpz_t b;
    mpz_t c;
    mpz_init(a);
    mpz_init(b);
    mpz_init(c);

    int N=48;
    int incr=4;
    int i=incr;
    int j;
    int reps=100;
    int k;

    for(;i<=N;i+=incr)
    {
        j=incr;
        for(;j<=i;j+=incr)
        {
            sprintf(name_in,"../data/data_in_%dx%d.dat",i,j);
            sprintf(name_out,"../data/data_out_%dx%d.dat",i,j);
            in=fopen(name_in,"w");
            out=fopen(name_out,"w");
            if(in==NULL)
            {
                printf("failed 2 open: %s",name_in);
            }
            if(out==NULL)
            {
                printf("failed 2 open: %s",name_out);
            }

            k=0;
            for(;k<reps;k++)
            {
                do //generate randomprime
                {
                    mpz_urandomb(a,rnd,i);
                }while(!mpz_probab_prime_p(a,128));
                gmp_fprintf(out,"%Zd\n",a);
                do
                {
                    mpz_urandomb(b,rnd,j);
                }while(!mpz_probab_prime_p(b,128));
                
                gmp_fprintf(out,"%Zd\n",b);
                
                gmp_fprintf(out,"\n");

                mpz_mul(c,a,b);
                gmp_fprintf(in,"%Zd\n",c);

            }
            printf("%dx%d\n",i,j);
            fclose(in);
            fclose(out);
        }
    }
    return 0;
    
}
