/* sieve.c - sieve part of the Multiple Polynomial Quadratic Sieve

  Copyright 1996 Scott Contini (see file README.orig)
  Copyright 2005 Paul Zimmermann (modified to use GMP).

  This program is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the
  Free Software Foundation; either version 2 of the License, or (at your
  option) any later version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
  more details.

  You should have received a copy of the GNU General Public License along
  with this program; see the file COPYING.  If not, write to the Free
  Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
  02111-1307, USA.
*/

#include "mpqs.h"

#define NumBits(a) mpz_sizeinbase (a,2)

#define BLOCK_LENGTH		60000	/* amount of sieve array to
					sieve upon at a time */
#define NUM_BLOCKS		10	/* M = BLOCK_LENGTH * NUM_BLOCKS */
#define START                   20	/* where to start sieving */
#define THRESH1                 50
#define ERROR2                  4

	/* note: we are assuming that any prime that divides the
	multiplier must have index less than START.  In other words,
	the primes that you sieve with are all larger than anything
	that divides the multiplier. */


char SMOOTH_FILE[1000], SEMI_FILE[1000], PP_FILE[1000], GOOD_PRIME_FILE[1000];
char AVAL_FILE[1000], STATUS_FILE[1000], ADATA_FILE[1000];


#define MAX_FB                  210000   /* Max primes in Factor Base */

#if (GMP_LIMB_BITS == 32)
   #define LONG long long
   #define LOG2_CHAR_PER_LONG 2
   #if (THRESH1 > 64)
   #define COMP                    0x80808080
   #else
   #define COMP                    0xc0c0c0c0
   #endif
#elif (GMP_LIMB_BITS == 64)
   #define LONG long
   #define LOG2_CHAR_PER_LONG 3
   #if (THRESH1 > 64)
   #define COMP                    0x8080808080808080
   #else
   #define COMP                    0xc0c0c0c0c0c0c0c0
   #endif
#else
#error "invalid value of GMP_BITS_PER_LIMB"
#endif

#define CHAR_PER_LONG (1 << LOG2_CHAR_PER_LONG)

int	comp;
int     rootN[MAX_FB];                  /* sqrt(N) mod primes[i] */
int	soln1[MAX_FB],soln2[MAX_FB];	/* for sieving */
int	ptr1[MAX_FB],ptr2[MAX_FB];	/* for sieving */
int     primes[MAX_FB],                 /* primes in factor base */
	p_count;
char	mprimes[START];			/* multiplier primes */
int	mp[10];				/* also multiplier primes */
int	mpc=0;
unsigned char   accum[BLOCK_LENGTH+1000];
unsigned char   log_p[MAX_FB];          /* logs of the primes in factor base */
mpz_t        N;			        /* number to factor */
mpz_t        Nby2;
mpz_t        a;
mpz_t        a2;
mpz_t        ainv;
mpz_t        b;
unsigned char   thresh2;
int     ssc,                            /* number of semi-smooth residues */
        smc,                            /* number of smooth residues */
        ppc;                            /* number of pps */
int     STOP = -1;                      /* When to stop sieving */
int     MACHINE_NO = 0;                 /* For parallel factoring */
unsigned long LP;                       /* large prime bound */
mpz_t semi_thresh;                      /* sieving threshold */
mpz_t pp_thresh;
mpz_t ztemp1, ztemp2, ztemp3;
int	FSM, FSS, FPP, F;		/* # smooths, semis, pps */
int	multiplier = 0;
mpz_t B2;
int ssc=0, smc=0, ppc=0;
int INIT;
int M;					/* length of sieve interval */
long factor_time = 0;                   /* time spent in large prime factor */

int smallp[14] = { 2,3,5,7,11,13,17,19,23,29,31,37,41,43};


FILE *fsmooth, *fsemi, *fadata, *fpp;

static int
mpz_jacobi_si_si (long j, long i)
{
  mpz_t J;
  int res;

  mpz_init_set_si (J, j);
  res = mpz_kronecker_si (J, i);
  mpz_clear (J);
  return res;
}

/* return factors of n (product of two large primes) into p[0] and p[1].
   If unsuccessful, return 0.
*/
static int
factor (unsigned long *p, mpz_t n)
{
  mpz_t x, y, t;
  unsigned long a = 1;
  int res = 0;

  mpz_init (x);
  mpz_init (y);
  mpz_init (t);
  mpz_set_ui (x, 1);
  mpz_set (y, x);
  while (1)
    {
      mpz_mul (x, x, x);
      mpz_add_ui (x, x, a);
      mpz_mod (x, x, n);
      mpz_mul (y, y, y);
      mpz_add_ui (y, y, a);
      mpz_mod (y, y, n);
      mpz_mul (y, y, y);
      mpz_add_ui (y, y, a);
      mpz_mod (y, y, n);
      mpz_sub (t, x, y);
      mpz_gcd (t, t, n);
      if (mpz_cmp_ui (t, 1))
	break;
    }
  res = mpz_cmp (t, n); /* non-zero if non trivial factor found */
  if (res && mpz_fits_ulong_p (t))
    {
      p[0] = mpz_get_ui (t);
      mpz_divexact_ui (t, n, p[0]);
      if (mpz_fits_ulong_p (t))
	p[1] = mpz_get_ui (t);
      else
	res = 0;
    }
  mpz_clear (x);
  mpz_clear (y);
  mpz_clear (t);
  return res;
}

static void
trial_divide (int index, int sd[100], int sdc)
{

	register int i,p;
	int e;
        mpz_t residue, u;
	int sign=1,j;
	static int all[100];
	int allc=0,mi;
	static int exp[100];
	int tsdc = sdc;

        mpz_init (residue);
        mpz_init (u);

	/* sd holds the index of the small prime divisors of the
	residue (excluding any primes that divide the multiplier).
	Now find the large primes that divide the residue. */

	for (i= START; i<p_count; i++) {
		register int test;

		p = primes[i];
		test = index % p;
		if (test < 0) test += p;
		if (test == soln1[i] || test == soln2[i])
			sd[tsdc++]=i;
	}

	/* compute the residue */

        mpz_mul_si (ztemp1, a2, index);
        mpz_add (ztemp2, ztemp1, b);
        mpz_mul (u, ztemp2, ainv);
        mpz_mod (u, u, N);
        mpz_mul (residue, u, u);
        mpz_mod (residue, residue, N);  /* probably can do this faster */
	if (mpz_cmp (residue, Nby2) > 0)
          {
            mpz_sub (residue, N, residue);
            sign = -1;
          }

	/* now trial divide and also check the primes which divide
	the multiplier.  Put all prime divisors into the all array.
	Also find their multiplicities and put that in exp. */

	for (i = 0, mi = 0; i < tsdc; ) {

		/* if there are still multiplier primes to check
		and the next is smaller than the next prime in the
		sd array... */
		if ((mi< mpc) && (mp[mi] < sd[i])) {
			p = primes[mp[mi]];
                        j = mpz_tdiv_q_ui (ztemp1, residue, p);
			if (j == 0) {
                          e = 1;
                          mpz_swap (residue, ztemp1); /* faster than mpz_set */
                          while (mpz_tdiv_q_ui (ztemp1, residue, p) == 0)
                            {
                              e++;
                              mpz_swap (residue, ztemp1); /* faster than set */
                            }
                          all[allc] = mp[mi]+1;
                          exp[allc++] = e;
			}

			mi++;
		}
		else {
			p = primes[sd[i]];
                        j = mpz_tdiv_q_ui (residue, residue, p);
			if (j != 0) {
				printf("ERROR:\n");
				printf("u=");
                                mpz_out_str (stdout, 10, u);
                                printf ("\n");
				printf("index=%d\n",index);
				printf("a=");
                                mpz_out_str (stdout, 10, a);
                                printf ("\n");
				printf("b=");
                                mpz_out_str (stdout, 10, b);
                                printf ("\n");
				printf("sieving says prime %d ",p);
				printf("should divide this residue\n");
				printf("remaining residue = ");
				mpz_mul_ui (residue, residue, p);
				mpz_add_ui (residue, residue, j);
                                mpz_out_str (stdout, 10, residue);
                                printf ("\n");
				printf("divisors so far:\n");
				for (j=0; j<allc; j++)
					printf("%d ",primes[all[j]]);
				printf("\n");
				exit(1);
			}
			e = 1;
			while (mpz_tdiv_q_ui (ztemp1, residue, p) == 0)
                          {
                            e++;
                            mpz_swap (residue, ztemp1); /* faster than set */
                          }
			all[allc] = sd[i]+1;
			exp[allc++] = e;
			i++;
		}
	}

	/* if the remaining residue is small enough, then keep it! */

	if (mpz_cmp_ui (residue, 1) == 0) {
		/* full relation */

		smc ++;
                mpz_out_str (fsmooth, 10, u);
		fprintf (fsmooth, " ");
		if (sign == -1)
                  fprintf (fsmooth, "-1 1 ");
		for (j = 0; j < allc; j++)
                  fprintf (fsmooth, "%d %d ", all[j], exp[j]);
		fprintf (fsmooth, "0\n");
	}
	else if (mpz_cmp_ui (residue, LP) < 0)
          {
		/* single large prime: necessarily prime since LP < B^2 */

		ssc ++;
                mpz_out_str (fsemi, 10, residue);
		fprintf (fsemi, " ");
                mpz_out_str (fsemi, 10, u);
		fprintf (fsemi, " ");
		if (sign == -1)
                  fprintf (fsemi, "-1 1 ");
		for (j = 0; j < allc; j++)
                  fprintf (fsemi, "%d %d ", all[j], exp[j]);
		fprintf (fsemi, "0\n");
	}
	else if (mpz_probab_prime_p (residue, 1) == 0)
          {
	    int res;
            /* partial-partial relation */
            unsigned long p[2];
            factor_time -= cputime ();
            res = factor (p, residue);
            factor_time += cputime ();
	    /* check both factors are less than large prime bound */
	    if (res != 0 && p[0] < LP && p[1] < LP)
	      {
		ppc ++;
		fprintf (fpp, "%lu %lu ", p[0], p[1]);
		mpz_out_str (fpp, 10, u);
		fprintf (fpp, " ");
		if (sign == -1)
		  fprintf (fpp, "-1 1 ");
		for (j = 0; j < allc; j++)
		  fprintf (fpp, "%d %d ", all[j], exp[j]);
		fprintf (fpp, "0\n");
	      }
          }

        mpz_clear (residue);
        mpz_clear (u);
}

static void
sieve (void)
{
  int k;
  int sdc, index;
  static int sd[100];		/* small divisors */

  /* first sieve to the right */

  for (k = 0; k < NUM_BLOCKS; k++)
    {
      register int i;
      register int p;
      register char logp;

      memset (accum, INIT, BLOCK_LENGTH + 1000);

      i = p_count - 1;
      p = primes[i];
      logp = log_p[i];

      for (;  i >= START ;  )
	{
	  register int i1, i2;

	  i1 = ptr1[i];
	  i2 = ptr2[i];

	  while (i2 < BLOCK_LENGTH)
	    {
	      accum[i1] += logp;
	      accum[i2] += logp;
	      i1 += p;
	      i2 += p;
	    }

	  if (i1 < BLOCK_LENGTH)
	    {
	      accum[i1] += logp;
	      i1 += p;
	    }

	  /* adjust ptrs */

	  ptr1[i] = i1 - BLOCK_LENGTH;
	  ptr2[i] = i2 - BLOCK_LENGTH;
	  if (ptr2[i] < ptr1[i])
	    {
	      p = ptr2[i];
	      ptr2[i] = ptr1[i];
	      ptr1[i] = p;
	    }

	  /* advance to next prime */
	  i--;
	  p = primes[i];
	  logp = log_p[i];
	}

      /* now scan the array for hits */

      {
	register long *fscan= (long*) accum;
	register long  *lim = (fscan + (BLOCK_LENGTH >> LOG2_CHAR_PER_LONG));

	for (; fscan <=  lim; fscan++ )
	  /* scan array CHAR_PER_LONG cells at a time */
	  if ((*fscan) & COMP) /* found hit */
	    {
	      register unsigned char *scan = (unsigned char*) fscan;
	      register int j;

	      for (j = 0; j < CHAR_PER_LONG; j++)
		if (scan[j] & comp)
		  {
		    int test;

		    sdc = 0;
		    index = (int) (scan-accum);
		    index = index + j + k * BLOCK_LENGTH;

		    for (i = 0; i < START; i++)
		      if (mprimes[i] == 0)
			{
			  p = primes[i];
			  test = index % p;
			  if (test == soln1[i] || test == soln2[i])
			    {
			      sd[sdc++] = i; /* save index */
			      scan[j] += log_p[i];
			    }
			}

		    if (scan[j] >= thresh2) /* try it */
		      {
			sd[sdc] = 0;
			trial_divide (index, sd, sdc);
		      }
		  }
	    }
      }
      
    }


	/* now sieve to the left */


	for (k=p_count-1; k>=0; k--) {
		/* we switch the order of the ptrs here so
		that ptr1 is less than ptr2 */
		ptr2[k] = primes[k] - soln1[k];
		ptr1[k] = primes[k] - soln2[k];
	}


	for (k=0; k < NUM_BLOCKS; k++) {
		register int i;
		register int p;
		register char logp;

		memset (accum, INIT, BLOCK_LENGTH+1000);

		i = p_count-1;
		p = primes[i];
		logp = log_p[i];

		for (;  i >= START ;  ) {
			register int i1,i2;

			i1 = ptr1[i];
			i2 = ptr2[i];

			while (i2 < BLOCK_LENGTH) {
				accum[i1] += logp;
				accum[i2] += logp;
				i1 += p;
				i2 += p;
			}

			if (i1 < BLOCK_LENGTH) {
				accum[i1] += logp;
				i1 += p;
			}

			/* adjust ptrs */

			ptr1[i] = i1 - BLOCK_LENGTH;
			ptr2[i] = i2 - BLOCK_LENGTH;
			if (ptr2[i] < ptr1[i]) {
				p = ptr2[i];
				ptr2[i] = ptr1[i];
				ptr1[i] = p;
			}

			/* advance to next prime */
			i--;
			p = primes[i];
			logp = log_p[i];
		}

		/* now scan the array for hits */

		{
		register long *fscan=( long *) accum;
		register long  *lim= (fscan+(BLOCK_LENGTH>>LOG2_CHAR_PER_LONG));

		for (; fscan <=  lim; fscan++ )
		  /* scan array CHAR_PER_LONG cells at a time */
		    if ((*fscan) & COMP) {		/* find hit */
			register unsigned char *scan =(unsigned char *)
			   fscan;
			register int j;

			for (j=0; j<CHAR_PER_LONG; j++)
			    if (scan[j] & comp) {
				int test;

				sdc=0;
				index= (int) (scan-accum);
				index= index + j + k*BLOCK_LENGTH;
				index= -index;

				for (i=0; i < START; i++)
				    if (mprimes[i]==0) {
					p = primes[i];
					test = index % p;
					if (test < 0) test += p;
					if (test== soln1[i] || test==
					 soln2[i]) {
					    sd[sdc++]= i; /* save index */
					    scan[j] +=log_p[i];
					}
				    }

				if (scan[j] >= thresh2) {	/* try it */
				    sd[sdc]=0;
				    trial_divide (index, sd, sdc);
				}
			   }
		    }
		}

	}


}

static int
inverse (int s, int t)
{
	register int u1,v1,q;
        register int u2,v2;
 
	u1=1; v1=0;
        u2=s; v2=t;

	while (v2 != 0) {
	  /* unroll twice and swap u/v */
		q = u2 / v2;
		u1 = u1 - q * v1;
		u2 = u2 - q * v2;

		if (u2 == 0)
		  {
		    u1 = v1;
		    break;
		  }

		q = v2 / u2;
		v1 = v1 - q * u1;
		v2 = v2 - q * u2;
        }
 
	if (u1 < 0)
		u1=u1-t*(-u1/t-1);
 
	return u1;
}

static void
compute_roots () {

	int i;
	int bmodp, amodp;

	for (i = 0 ; i < p_count; i++) {
		register int p;
		register LONG x, inv;

		if (i < START && mprimes[i]==1)
			continue;

		p = primes[i];
                amodp = mpz_tdiv_ui (a, p);
                bmodp = mpz_tdiv_ui (b, p);

		x = rootN[i] - bmodp;
		if (x < 0) x += p;
		inv = inverse (amodp, p);
		soln1[i] = (((x * inv) % p) * inv) % p;

		x = p - rootN[i] - bmodp;
		if (x < 0) x += p;
		soln2[i] = (((x * inv) % p) * inv) % p;

		if (soln2[i] < soln1[i]) {	/* swap 'em */
			p = soln1[i];
			soln1[i] = soln2[i];
			soln2[i] = p;
		}
		ptr1[i] = soln1[i];
		ptr2[i] = soln2[i];

	}
}

/* Function: expmod
 
   Description: calculates and returns a^2 mod n.  Use the simple,
        recursive algorithm.
*/
static int
expmod (LONG a, LONG e, LONG n)
{
  LONG temp;
 
  if (e == 1)
    return a % n;

  if (e & 1)
    return (a * expmod (a, e-1, n)) % n;

  temp= expmod (a, e >> 1, n);
  return (temp * temp) % n;
}

/* Function: jacobi
 
   Description: calculates and returns the jacobi symbol.
 
   Note: The jacobi symbol is the same as the Legendre symbol mod a prime
*/
static int
jacobi (int m, int n)
{
 
	int t;
 
	if (m>n) m=m%n;
	if (m==0) return 0;
	if (m==1) return 1;
	if (m==2) {
		if ((n & 1) == 0) return 0;
		t=n&7;                                  /* t = n % 8 */
		if (t==1 || t==7) return 1;
		return -1;
	}
	if ((m&1)==0)
		return jacobi((m>>1),n)*jacobi(2,n);
	if ((m&3)==3 && (n&3)==3) return -jacobi(n,m);
	return jacobi(n,m);
	 
}

/* Function: sqrtmod
 
   Description: Calculates the squareroot of a mod prime p.
	Uses the RESSOL algorithm of Shanks.  For description of algorithm,
	see:  An Introduction to the Theory of Numbers, 5th edition.
	Authors: Niven, Zuckerman, and Montgomery.  Pages 110-114.
 
   Always return the smallest squareroot
*/
static int
sqrtmod (int a, int p)
{
 
	register LONG r,s,c;
	LONG m,n,k,kprime,i,j;
 
	/* if number is small then search for solution */
 
	if ((p < 50) || ((p < 200) && ((p & 2) == 0)))
		for (i=1; i<p; i++)
			if (i*i % p == a) return i;
 
	if ((p & 2) == 1) {   /* p == 3 mod 4 */
		s=expmod((LONG) a,(LONG) ((p+1)>>2),(LONG) p);
		if (s > (p>>1)) return p-s;
		else return s;
	}
 
	kprime= -1; k=0;
	m=p-1;
	while ((m&1)==0) k++, m>>=1;
 
	/* perform two fast exponentiations at the same time below */
 
	s=(m+1)>>1; c=a; r=1; n=1;
	for (i=0; (j=(1<<i)) <= m; i++) {
		if (m&j) n=(n*c) % p;
		if (s&j) r=(r*c) % p;
		c=(c*c)%p;
        }
 
	while(1) {
		if (n == 1)
                  {
                    if (r > (p>>1)) return p-r;
                    else return r;
                  }
		if (kprime < 0) {
			for (i=0; ;i++) {
				c=primes[i];
				if (jacobi((int) c,p)==-1) break;
			}
			c=expmod(c,m,(LONG) p);
		}
		s=n;
		for (kprime=1; ; kprime++) {
			s=(s*s)%p;
			if (s==1) break;
		}
		s=c;
		for (i=0; i<k-kprime-1; i++)
			s=(s*s) % p;
		r=(s*r) % p; c=(s*s) % p; n=(c*n) % p; k=kprime;
	}
 
 
}

/******************************************************************
 * computes x so that x^2 == a mod p for prime p, and puts x in s.
 * if no such x exists or if p is not prime, then sets s=0.             
 ******************************************************************/
static void
mpz_sqrtmod (mpz_t s, mpz_t a, mpz_t p)
{

	mpz_t q, nn, x, y, t, b, z;
	register long n, i, m, r;

        mpz_init (q);
        mpz_init (y);

	if (mpz_fdiv_ui (p, 4) == 3)
          {
            mpz_add_ui (q, p, 1);
            mpz_div_2exp (q, q, 2);
            mpz_powm (s, a, q, p);
            mpz_mul (y, s, s);
            mpz_mod (y, y, p);
            if (mpz_cmp (y, a))
              mpz_set_ui (s, 0);
            mpz_clear (q);
            mpz_clear (y);
            return;
          }

        mpz_init (nn);
        mpz_init (x);
        mpz_init (t);
        mpz_init (b);
        mpz_init (z);

        mpz_sub_ui (q, p, 1);
        if (mpz_cmp_ui (q, 0) == 0)
          r = -1;
        else
          {
            r = mpz_scan1 (q, 0);
            mpz_div_2exp (q, q, r);
          }
	n = 3;
	if ((mpz_fdiv_ui (p, 4) & n) == 3)
          m = -mpz_kronecker_si (p, n);
	else
          m = mpz_kronecker_si (p, n);

	while (m != -1) {
loop:
		n += 2;
		if (!(n % 3))
			n += 2;
		for (i = 5; i * i <= n; )  {
			if (!(n % i))
                          goto loop;
			i += 2;
			if (!(n % i))
                          goto loop;
			i += 4;
		}
		if ((mpz_fdiv_ui (p, 4) & n) == 3)
                  m = -mpz_kronecker_si (p, n);
		else
                  m = mpz_kronecker_si (p, n);
	}

        mpz_set_si (nn, n);
        
        mpz_powm (z, nn, q, p);

        mpz_set (y, z);
        mpz_sub_ui (t, q, 1);
        mpz_div_2exp (t, t, 1); /* t = (q-1)/2 */
        mpz_powm (x, a, t, p);
        mpz_mul (t, x, x);
        mpz_mod (t, t, p);
        mpz_mul (b, t, a);
        mpz_mod (b, b, p);      /* b = a*x^2 mod p */
        mpz_mul (x, a, x);
        mpz_mod (x, x, p);

	while (1) {
		if (mpz_cmp_ui (b, 1) == 0)
                  {
                    mpz_set (s, x);
                    goto end;
                  }

                mpz_mul (t, b, b);
                mpz_mod (t, t, p);
		m = 1;
		while (mpz_cmp_ui (t, 1) != 0) 
                  {
                    if (m == r)
                      {
                        mpz_set_ui (s, 0);
                        goto end;
                      }
                    mpz_mul (t, t, t);
                    mpz_mod (t, t, p);
                    m ++;
                  }

                mpz_set (t, y);
		if (r <= m)
                  {
                    mpz_set_ui (s, 0);
                    goto end;
                  }
		for (i = r - m - 1; i; i--)
                  {
                    mpz_mul (t, t, t);
                    mpz_mod (t, t, p);
                  }
                mpz_mul (y, t, t);
                mpz_mod (y, y, p);
		r = m;
                mpz_mul (x, x, t);
                mpz_mod (x, x, p);
                mpz_mul (b, b, y);
                mpz_mod (b, b, p);
	}
 end:
        mpz_clear (q);
        mpz_clear (nn);
        mpz_clear (x);
        mpz_clear (y);
        mpz_clear (t);
        mpz_clear (b);
        mpz_clear (z);
}

/*
	Find res so that   res * res == x mod p^2, where
	p is prime.  The input asks for p2 which should be p^2.
*/
static void
zsqrtp2 (mpz_t x, mpz_t p, mpz_t p2, mpz_t res)
{
  mpz_t rr, xmodp, xmodp2; /* FIXME: use static variables? */

  mpz_init (rr);
  mpz_init (xmodp);
  mpz_init (xmodp2);

  mpz_mod (xmodp, x, p);
  mpz_mod (xmodp2, x, p2);

  mpz_sqrtmod (rr, xmodp, p);              /* rr^2 == x mod p */

  mpz_mul (ztemp1, rr, rr);
  mpz_sub (ztemp1, ztemp1, xmodp2); /* ztemp1 is rr^2 - x mod p^2 */
  mpz_mul_ui (ztemp3, rr, 2);
  mpz_invert (ztemp2, ztemp3, p);   /* ztemp2 is (2*ztemp1)^-1 mod p */
  mpz_mul (ztemp3, ztemp1, ztemp2);
  mpz_sub (ztemp1, rr, ztemp3);
  mpz_mod (res, ztemp1, p2);

  mpz_clear (rr);
  mpz_clear (xmodp);
  mpz_clear (xmodp2);
}

/*
	This routine just searches for the next prime after
	a which has the property that N is a quadratic residue
	mod it.  Once that is found, this becomes the new value of
	a and b is computed so that b^2 == N mod a^2, and b is between
	0 and a^2/2.
*/
static void
get_coeffs ()
{

	static char array[99];
	register int i,j,k,p;
	FILE *fp;

loop:

	/* sieve for the next prime */
	memset(array, 0, 99);

        mpz_add_ui (a, a, 1);   /* start sieving at one more than
				previous a */

	/* trial divide by smallp array */
	for (k=0; k < 14; k++) {
		p = smallp[k];
                j = mpz_tdiv_ui (a, p); /* j is remainder after
                                           division by p */
		if (j==0) array[0] = 1;
		for (i= p-j; i<99; i+= p)
			array[i] = 1;
	}

	for (i=0; i < 99; i++) {
		if (array[i]) continue;
                mpz_add_ui (ztemp1, a, i);
		if (mpz_probab_prime_p (ztemp1, 1))
                  {
                    if (mpz_jacobi (N, ztemp1) == 1)
                      break;
                  }
	}

	if (i == 99)
          {			/* nothing found */
            mpz_add_ui (a, a, 99);
            goto loop;		/* try again */
          }

	/* now we have a prime for which N is a q.r. mod */
	/* get b: the sqrt of N mod a^2 */

        mpz_set (a, ztemp1);
        mpz_mul (a2, a, a);             /* a2 = a * a */
        mpz_invert (ainv, a, N);        /* a^-1 mod N */
	zsqrtp2 (N, a, a2, b);

	/* make sure b is between 0 and a^2/2 */

        mpz_sub (ztemp1, a2, b);
	if (mpz_cmp (ztemp1, b) < 0)
          mpz_set (b, ztemp1);

	/* record this a */

	fp = fopen (AVAL_FILE, "w");
        mpz_out_str (fp, 10, a);
	fprintf (fp,"\n");
	fclose (fp);

        mpz_out_str (fadata, 10, a);
	fprintf (fadata, "\n");

}

static void
usage ()
{
  fprintf (stderr, "Usage: sieve -p <params> -n <nfulls>\n");
  fprintf (stderr, "-m <machine_number>: for parallel sieving\n");
  fprintf (stderr, "-k <multiplier>: fix multiplier\n");
  fprintf (stderr, "-T <value>: use B^T as large prime bound [default 2.5]\n");
  exit (1);
}

int
main (int argc, char *argv[])
{
  FILE *fp, *status_fp;
  unsigned long i;
  double dthresh;
  int loopc=0;
  char *params = NULL;
  unsigned long B; /* largest prime in factor base */
  double LPn; /* number of large primes */
  double S1, S2;
  double T = 2.5;

  while (argc > 1)
    {
      if (argc > 2 && strcmp (argv[1], "-p") == 0)
        {
          params = argv[2];
          argv += 2;
          argc -= 2;
        }
      else if (argc > 2 && strcmp (argv[1], "-k") == 0)
        {
          multiplier = atoi (argv[2]);
          argv += 2;
          argc -= 2;
        }
      else if (argc > 2 && strcmp (argv[1], "-m") == 0)
        {
          MACHINE_NO = atoi (argv[2]);
          argv += 2;
          argc -= 2;
        }
      else if (argc > 2 && strcmp (argv[1], "-n") == 0)
        {
          STOP = atoi (argv[2]);
          argv += 2;
          argc -= 2;
        }
      else if (argc > 2 && strcmp (argv[1], "-T") == 0)
        {
          T = atof (argv[2]);
          argv += 2;
          argc -= 2;
        }
      else
        usage ();
    }

  if (params == NULL || STOP == -1)
    usage ();

#define GMP_VERSION __GNU_MP_VERSION ## . ## __GNU_MP_VERSION_MINOR ## . ## __GNU_MP_VERSION_PATCHLEVEL
  /* check GMP library vs headers */
  {
    char buffer[100];
    if (__GNU_MP_VERSION_PATCHLEVEL != 0)
      sprintf (buffer, "%d.%d.%d", __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR,
	       __GNU_MP_VERSION_PATCHLEVEL);
    else
      sprintf (buffer, "%d.%d", __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR);
    if (strcmp (gmp_version, buffer))
      {
	fprintf (stderr, "GMP header and library do not match: %s %s\n",
		 gmp_version, buffer);
	exit (1);
      }
  }

        mpz_init (N);
        mpz_init (Nby2);
        mpz_init (ztemp1);
        mpz_init (ztemp2);
        mpz_init (ztemp3);
        mpz_init (semi_thresh);
        mpz_init (B2);
        mpz_init (pp_thresh);
        mpz_init (a);
        mpz_init (a2);
        mpz_init (ainv);
        mpz_init (b);

	printf ("Welcome to Multiple Polynomial Quadratic Sieve\n");
	printf ("using THRESH1=%d ERROR2=%d\n", THRESH1, ERROR2);

	/* read number to factor from params file */

	read_params (N, &B, &LP, params);
	printf ("trying to factor:\n  ");
        mpz_out_str (stdout, 10, N);
        printf ("\n");

	M = BLOCK_LENGTH * NUM_BLOCKS;
	printf ("using M=%d (%d blocks of length %d)\n", M, NUM_BLOCKS,
                BLOCK_LENGTH);
	printf ("factor base bound %lu\n", B);
        printf ("large prime bound %lu\n", LP);

        /* sieving threshold is B^T */
        {
          double t;
          t = pow ((double) B, T);
          mpz_set_d (semi_thresh, t);
        }
        printf ("sieve threshold is ");
        mpz_out_str (stdout, 10, semi_thresh);
        printf ("\n");

        LPn = (double) LP;
        S1 = log (log (LPn)) - log (log ((double) B));
        S2 = 1.0 / (double) B / log ((double) B) - 1.0 / LPn / log (LPn);
        S2 = S2 / S1 / S1 / 2.0;
        LPn = LPn / log (LPn); /* estimation of # of primes < LP */
        LPn -= (double) B / log ((double) B); /* subtract primes < B */
        LPn = LPn / 2.0;       /* only about half the primes may appear */

        if (MACHINE_NO)
          printf ("Machine number: %d\n", MACHINE_NO);
	printf ("stop sieving after %d smooths\n", STOP);

        /* factor base is stored in <params>.fbase */
        strcpy (GOOD_PRIME_FILE, params);
        strcat (GOOD_PRIME_FILE, FBASE);
	printf("factor base primes will be put in file '%s'\n",
	   GOOD_PRIME_FILE);

        /* fulls are stored in <params>.fulls */
        strcpy (SMOOTH_FILE, params);
        strcat (SMOOTH_FILE, FULLS);
	printf("smooths will be put in file '%s'\n", SMOOTH_FILE);

        /* partials are stored in <params>.partials */
        strcpy (SEMI_FILE, params);
        strcat (SEMI_FILE, PARTIALS);
	printf("partials will be put in file '%s'\n", SEMI_FILE);

        /* 2-large prime partials are stored in <params>.pps */
        strcpy (PP_FILE, params);
        strcat (PP_FILE, PPS);
	printf ("pps will be put in file '%s'\n", PP_FILE);

        strcpy (AVAL_FILE, params);
        strcat (AVAL_FILE, AVAL);

        strcpy (STATUS_FILE, params);
        strcat (STATUS_FILE, STATUS);

        strcpy (ADATA_FILE, params);
        strcat (ADATA_FILE, ADATA);

	/* open output data files */

	fsmooth = fopen (SMOOTH_FILE, "a");
	fsemi = fopen (SEMI_FILE, "a");
	fpp = fopen (PP_FILE, "a");
	fadata = fopen (ADATA_FILE, "a");

        /* find best multiplier */
        if (multiplier == 0)
          multiplier = find_multiplier (N, (double) B);
        printf ("using multiplier k = %d\n", multiplier);
        mpz_mul_ui (N, N, multiplier);
        printf ("k*N = ");
        mpz_out_str (stdout, 10, N);
        printf ("\n");

	/* compute half of N */

        mpz_div_2exp (Nby2, N, 1);

	/* compute sieving thresholds */

        mpz_mul (B2, ztemp1, ztemp1);
        mpz_mul_ui (pp_thresh, B2, 64 * 64); /* largest partial partials that
                                                we will try to factor */
	

	/* The following code calculates the sieve threshold for the given
	parameters.  It then adjust the sieve threshold so that it is a power
	of 2.  This adjustment allows us to scan the sieve array much faster
	when we search for probable smooth residues. */

	dthresh = 0.5 * NumBits (N) * log (2.0) + log ((double) M)
          - NumBits (semi_thresh) * log (2.0) - ERROR2;
	thresh2 = (unsigned char) dthresh;       /* sieving threshold */
	printf ("sieve thresholds: %d %d\n", THRESH1, thresh2);
	fflush (stdout);
	INIT = 128;
	while (!(INIT&THRESH1)) INIT>>=1;
	if (INIT != THRESH1) INIT<<=1;	/* INIT>=THRESH, INIT is power of 2 */
	INIT = INIT - THRESH1;
	comp=0;
	for (i=THRESH1+INIT; i<=128; i<<=1)
		comp += i;
	if (comp != (COMP& 0xff)) {
		printf("please recompile with COMP = 0x%x%x%x%x\n",
			comp,comp,comp,comp);
		exit(1);
	}
	thresh2 += INIT;
	printf("modified: %d %d %d\n",INIT,THRESH1+INIT,thresh2);

	/* a must be size sqrt (sqrt(2*N)/M).  The first time the
	program is ran, we need to compute this, and store it in
	the a variable (with some adjustments to it if we parallelize).
	a represents the previous leading coefficient.  We are pretending
	that the previous leading coefficient was sqrt(sqrt(2*N)/M).  The
	next value of a will be the next prime after this.
	If the program has already been ran, then just read in a from
	the "a_val" file. */

	fp = fopen (AVAL_FILE, "r");
	if (fp == NULL) {
                mpz_mul_ui (ztemp1, N, 2);
                mpz_sqrtrem (ztemp2, ztemp3, ztemp1); /* ztemp2 = sqrt(2*N) */
                mpz_div_ui (ztemp1, ztemp2, M); /* ztemp1 = sqrt(2*N)/M */
                mpz_sqrtrem (a, ztemp2, ztemp1); /* a = sqrt(sqrt(2*N)/M) */
                mpz_sqrtrem (ztemp1, ztemp2, a);
                mpz_div_ui (ztemp1, ztemp1, 100); /* ztemp1 will be the adjustment
                                                     "distance" for parallelizing. */
                mpz_mul_ui (ztemp2, ztemp1, MACHINE_NO);
                mpz_add (a, a, ztemp2); /* this should prevent overlap */
	}
	else
          {
            mpz_inp_str (a, fp, 10);
            fclose (fp);
          }

	/* read in the factor base primes */

	fp = fopen (GOOD_PRIME_FILE, "r");
	if (fp != NULL) {
		int j;
		printf("reading in factor base primes!\n"); fflush(stdout);
		p_count=0;
		while(fscanf(fp, "%d", &primes[p_count]) > 0) {
			log_p[p_count]=(unsigned char) (log((double)
			   primes[p_count]) +.5);
			i = primes[p_count];
			if (i == 2) {
				if (multiplier %2 ==1) {
					rootN[0]=1;
					mprimes[0] = 0;
				}
				else {
					rootN[0]=0;
					mprimes[0]= 1;
					mp[mpc++] = 0;
				}
				p_count++;
				continue;
			}
			/* if prime divides the multiplier then we
			must treat it a different way */
			j = multiplier % i;
			if (j==0 && p_count >= START) {
				printf("Error: multiplier too large.\n");
				exit(1);
			}
			if (j==0) {
				mprimes[p_count]=1;
				mp[mpc++]=p_count;
			}
			else {
				if (p_count < START)
					mprimes[p_count] = 0;
				/* calculate the sqrt of N mod each prime in
				the factor base */
                                j = mpz_tdiv_ui (N, i);
				if (j == 0) {
                                  printf("N is divisible by %d\n",j);
                                  exit(1);
				}
				rootN[p_count]=sqrtmod(j, i);
			}
			p_count++;
		}
	printf("done reading factor base primes\n"); fflush(stdout);
	fclose(fp);
	}
	else {
		printf ("creating factor base primes file!\n");
		fp = fopen (GOOD_PRIME_FILE, "w");
		primes[0] = 2;
		fprintf (fp,"2\n");
		log_p[0] = (unsigned char) (log((double) 2.0) +.5);
		if (multiplier % 2 == 1) {
			rootN[0]=1;
			mprimes[0] = 0;
		}
		else {
			rootN[0]=0;
			mprimes[0]= 1;
			mp[mpc++]=0;
		}
		p_count=1;
                getprime (0.0);
		while ((i = (int) getprime (1.0)) < B )
                  {
                    register int j, js;
                    if (i == 2)
                      {
                        fprintf (stderr, "need larger PRIM_BND for this factor base\n");
                        exit (1);
                      }
                    j = mpz_tdiv_ui (N, i);
                    if (j == 0)
                      {
                        if ((multiplier % i) != 0)
                          {
                            fprintf (stderr, "N is divisible by %lu\n", i);
                            exit (1);
                          }
                        if (p_count >= START)
                          {
                            fprintf (stderr, "Error: multiplier too large.\n");
                            exit (1);
                          }
                        continue; /* don't use the multiplier */
                      }
                    else if ( (js = mpz_jacobi_si_si (j, i)) == -1)
                      continue;
                    primes[p_count] = i;
                    fprintf (fp, "%lu\n", i);
                    fflush (fp);
                    if (p_count < START)
                      {
                        if (j == 0)
                          {
                            mp[mpc++] = p_count;
                            mprimes[p_count] = 1;
                          }
                        else
                          {
                            mprimes[p_count] = 0;
                          }
                      }
                    rootN[p_count] = sqrtmod (j, i);
                    log_p[p_count] = (unsigned char) (log ((double) i) + 0.5);
                    p_count++;
                  }
		fclose (fp);
	}

	primes[p_count]=0;

	/* get data from status file */

	status_fp = fopen (STATUS_FILE, "r");
	if (status_fp != NULL) {
		fscanf(status_fp,"%d",&F);              /* # good primes */
		fscanf(status_fp,"%d",&FPP);            /* # pps */
		fscanf(status_fp,"%d",&FSS);            /* # semi */
		fscanf(status_fp,"%d",&FSM);            /* # smooth */
		fclose(status_fp);
        }
	else {
		status_fp = fopen (STATUS_FILE, "w");
		fprintf(status_fp,"%d 0 0 0\n",p_count);
		F=p_count;
		FPP=0;
		FSS=0;
		FSM=0;
		fclose(status_fp);
        }

	printf ("sieving begins at %d\n",primes[START]);
        fflush (stdout);

	system ("date");
        system ("hostname");


	/* start main loop */

	while (FSM < STOP) {
		ppc = 0;
		ssc = 0;
		smc = 0;

                /* from page 42 of Emmanuel Thome' PhD thesis:
                   if the large primes have probability p_i,
                   then the expected number of cycles for n partials
                   is n^2*S2/2, where S2 = sum(p_i^2).
                   
                   If we take the large primes between B and L,
                   let S1 = sum(1/p, p=B..L), then p_i = 1/p/S1,
                   thus S2 = sum(1/p^2, p=B..L)/S1
                   We have S1 ~ log(log(L)) - log(log(B))
                   and sum(1/p^2, p=B..L) ~ 1/B/log(B) - 1/L/log(L) */

		if ((loopc & 63) == 0)
                  {
                    double d = (double) FSS;
                    d = d * d * S2;
                    printf ("time=%1.1fs: %d ffs, %d fps, %d pps, est. %d cycles (%d wanted)\n",
                            (double) cputime () / 1e3, FSM, FSS, FPP, (int) d, F);
                    fflush (stdout);
                  }

		get_coeffs ();		/* get new value of a, b */

		compute_roots();

		sieve ();

		/* update smooth counts, etc... */

		FPP += ppc;
		FSS += ssc;
		FSM += smc;
		status_fp = fopen (STATUS_FILE, "w");
		fprintf(status_fp,"%d %d %d %d\n",F,FPP,FSS,FSM);
		fclose(status_fp);

		fflush (fsemi);
		fflush (fsmooth);
		fflush (fpp);

		loopc++;
	}

        fprintf (stderr, "sieving took %1.1fs (factor %1.1fs)\n",
		 (double) cputime () / 1e3, (double) factor_time / 1e3);
	fprintf (stderr, "%d smooth equations\n%d semi-smooth equations\n",
                 FSM, FSS);
	fprintf (stderr, "%d partial partials\n", FPP);
        
	fclose (fsmooth);
	fclose (fsemi);
	fclose (fpp);
	fclose (fadata);

        mpz_clear (N);
        mpz_clear (Nby2);
        mpz_clear (ztemp1);
        mpz_clear (ztemp2);
        mpz_clear (ztemp3);
        mpz_clear (semi_thresh);
        mpz_clear (B2);
        mpz_clear (pp_thresh);
        mpz_clear (a);
        mpz_clear (a2);
        mpz_clear (ainv);
        mpz_clear (b);

        getprime (0.0); /* free prime buffer */

        return 0;
}
