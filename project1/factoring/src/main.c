#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#define NUM_INPUT 100

#define MAX_FACTORS 512

#define MAX_NUMBER_OF_TRIES 260000 //190000 (with 300 (and with 900) trial division) //202000 // < 205 000

#define NEXT_VALUE(r, x, n, c) \
    mpz_set_ui(cp, c);          \
    mpz_addmul(cp, x, x);       \
    mpz_mod(r, cp, n)

#define IS_ONE(x)  mpz_cmp_ui(x, 1)
#define IS_PRIME(x) mpz_probab_prime_p(x, 3)
#define FAILED(x) mpz_cmp_ui(x, 0) == 0
#define ADD_FACTOR(factors, num_factors, factor) \
    do {                                         \
        mpz_set( factors[num_factors], factor);  \
        num_factors++;                           \
    } while (0)

#define RAND ((int) (((double) rand()) / ((double)RAND_MAX)) * 100)
#include "primes.h"

mpz_t cp;

#include "primes.h"

void pollard_brent(mpz_t x_1, mpz_t x_2, mpz_t d, mpz_t n, int init, int c, int counter)
{
    int i;
    for(i = 0; i < 10000; ++i) {
        if(mpz_divisible_ui_p(n, primes[i])) {
            mpz_set_ui(d, primes[i]);
            return;
        }
    }

    unsigned long power = 1;
    unsigned long lam = 1;

    mpz_set_ui(x_1, init);

    mpz_set_ui(d, 1);
    while(mpz_cmp_ui(d, 1) == 0) {
        if (power == lam) {
            mpz_set(x_2, x_1);
            power *= 2;
            lam = 0;
        }

        mpz_mul(x_1, x_1, x_1);
        mpz_add_ui(x_1, x_1, c);
        mpz_mod(x_1, x_1, n);

        mpz_sub(d, x_2, x_1);
        mpz_abs(d, d);
        mpz_gcd(d, d, n);

        lam++;
        
        if (counter < 0){
            mpz_set_ui(d, 0);
            return;
        } else {
            --counter;
        }
    }

    if(mpz_cmp(d, n) == 0) {
        return pollard_brent(x_1, x_2, d, n, rand(), rand(), counter);
    }
}

void pollard(mpz_t x_1, mpz_t x_2, mpz_t d, mpz_t n, int init, int c, int counter)
{
    if(mpz_divisible_ui_p(n, 2)) {
        mpz_set_ui(d, 2);
        return;
    }

    mpz_set_ui(x_1, init);
    NEXT_VALUE(x_1, x_1, n, c);
    NEXT_VALUE(x_2, x_1, n, c);

    mpz_set_ui(d, 1);
    while (mpz_cmp_ui(d, 1) == 0) {
        NEXT_VALUE(x_1, x_1, n, c);
        NEXT_VALUE(x_2, x_2, n, c);
        NEXT_VALUE(x_2, x_2, n, c);

        mpz_sub(d, x_2, x_1);
        mpz_abs(d, d);
        mpz_gcd(d, d, n);
    
        if (counter < 0){
            mpz_set_ui(d, 0);
            return;
        } else {
            --counter;
        }
    }

    if(mpz_cmp(d, n) == 0) {
        return pollard(x_1, x_2, d, n, primes[RAND], primes[RAND], counter);
    }
}

int factor(mpz_t x_1, mpz_t x_2, mpz_t d, mpz_t n, mpz_t tmp, mpz_t *factors)
{
    int num_factors = 0;

    while (IS_ONE(n) != 0 && IS_PRIME(n) == 0) {
        /*pollard(x_1, x_2, d, n, 2, 1, MAX_NUMBER_OF_TRIES);*/
        pollard_brent(x_1, x_2, d, n, 2, 3, MAX_NUMBER_OF_TRIES);

        if(FAILED(d)) {
            return 0;
        }

        while (IS_ONE(d) != 0 && IS_PRIME(d) == 0) {
            /*pollard(x_2, x_2, tmp, d, 2, 1, MAX_NUMBER_OF_TRIES);*/
            pollard_brent(x_1, x_2, d, n, 2, 3, MAX_NUMBER_OF_TRIES);

            if(FAILED(tmp)) {
                return 0;
            }

            mpz_set(d, tmp);
        }
        do {
            ADD_FACTOR(factors, num_factors, d);
            mpz_divexact(n, n, d);
        } while (mpz_divisible_p(n, d));
    }
    if(IS_ONE(n) != 0) {
        ADD_FACTOR(factors, num_factors, n);
    }
    return num_factors;
}

void print_factors(mpz_t *factors, int num_factors)
{
    int i;

    if (num_factors == 0) {
        printf("fail\n");
    } else {
        for (i = 0; i < num_factors; ++i) {
            gmp_printf("%Zd\n", factors[i]);
        }
    }
    putchar('\n');
}

int main(void)
{
    int i, num_factors;
    mpz_t n, d, tmp, x_1, x_2;
    mpz_init(n); 
    mpz_init(d);
    mpz_init(tmp);
    mpz_init(x_1);
    mpz_init(x_2);
    mpz_init(cp);

    srand(1337);

    mpz_t *factors = malloc(sizeof(mpz_t) * MAX_FACTORS);

    for(i = 0; i < MAX_FACTORS; ++i) {
        mpz_init(factors[i]);
    }

    for (i = 0; i < NUM_INPUT; ++i) {
        mpz_inp_str(n, stdin, 10);
        num_factors = factor(x_1, x_2, d, n, tmp, factors);
        print_factors(factors, num_factors);
    }

    return 0;
}
