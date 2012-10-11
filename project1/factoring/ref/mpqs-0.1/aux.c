/* 
  Copyright 2005 Paul Zimmermann.

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

/* copied from ecm-6.0.1, file getprime.c */
double
getprime (double pp)
{
  static double offset = 0.0; /* offset for current primes */
  static long int current = -1; /* index of previous prime */
  static unsigned *primes = NULL; /* table of small primes up to sqrt(p) */
  static unsigned long int nprimes = 0; /* length of primes[] */
  static unsigned char *sieve = NULL; /* sieving table */
  static long int len = 0; /* length of sieving table */
  static unsigned long int *moduli = NULL; /* offset for small primes */

  if (pp == 0.0) /* free the tables, and reinitialize */
    {
      offset = 0.0;
      current = -1;
      free (primes);
      primes = NULL;
      nprimes = 0;
      free (sieve);
      sieve = NULL;
      len = 0;
      free (moduli);
      moduli = NULL;
      return pp;
    }

  while ((++current < len) && (sieve[current] == 0));

  if (current < len) /* most calls will end here */
    return offset + 2.0 * (double) current;

  /* otherwise we have to sieve */
  offset += 2.0 * (double) len;

  /* first enlarge sieving table if too small */
  if ((double) len * (double) len < offset)
    {
      free (sieve);
      len *= 2;
      sieve = (unsigned char *) malloc (len * sizeof (unsigned char));
      /* assume this "small" malloc will not fail in normal usage */
    }

  /* now enlarge small prime table if too small */
  if ((nprimes == 0) || (primes[nprimes-1] < sqrt(offset + len)))
      {
	if (nprimes == 0) /* initialization */
	  {
	    nprimes = 1;
	    primes = (unsigned *) malloc (nprimes * sizeof(unsigned long int));
	    /* assume this "small" malloc will not fail in normal usage */
	    moduli = (long unsigned int *) malloc (nprimes *
                                                   sizeof(unsigned long int));
	    /* assume this "small" malloc will not fail in normal usage */
	    len = 1;
	    sieve = (unsigned char *) malloc(len *
                                       sizeof(unsigned char)); /* len=1 here */
	    /* assume this "small" malloc will not fail in normal usage */
	    offset = 5.0;
	    sieve[0] = 1; /* corresponding to 5 */
	    primes[0] = 3;
	    moduli[0] = 1; /* next odd multiple of 3 is 7, i.e. next to 5 */
	    current = -1;
	    return 3.0;
	  }
	else
	  {
	    unsigned long int i, p, j, ok;

	    i = nprimes;
	    nprimes *= 2;
	    primes = (unsigned *) realloc (primes, nprimes *
                                           sizeof(unsigned long int));
	    moduli = (unsigned long int *) realloc (moduli, nprimes *
                                                    sizeof(unsigned long int));
	    /* assume those "small" realloc's will not fail in normal usage */
	    for (p = primes[i-1]; i < nprimes; i++)
	      {
		/* find next (odd) prime > p */
		do
		  {
		    for (p += 2, ok = 1, j = 0; (ok != 0) && (j < i); j++)
		      ok = p % primes[j];
		  }
		while (ok == 0);
		primes[i] = p;
		/* moduli[i] is the smallest m such that offset + 2*m = k*p */
		j = (unsigned long) fmod (offset, (double) p);
		j = (j == 0) ? j : p - j; /* -offset mod p */
		if ((j % 2) != 0)
		  j += p; /* ensure j is even */
		moduli[i] = j / 2;
	      }
	  }
      }

  /* now sieve for new primes */
  {
    long int i;
    unsigned long int j, p;
    
    for (i = 0; i < len; i++)
      sieve[i] = 1;
    for (j = 0; j < nprimes; j++)
      {
	p = primes[j];
	for (i = moduli[j]; i < len; i += p)
	  sieve[i] = 0;
	moduli[j] = i - len; /* for next sieving array */
      }
  }

  current = -1;
  while ((++current < len) && (sieve[current] == 0));

  return offset + 2.0 * (double) current;
}

unsigned long
count_lines (char *f)
{
  FILE *fp;
  int c;
  unsigned long n = 0;

  fp = fopen (f, "r");
  if (fp == NULL)
    {
      fprintf (stderr, "Error, cannot open %s\n", f);
      exit (1);
    }
  while ((c = fgetc (fp)) != EOF)
    n += c == '\n';
  fclose (fp);
  return n;
}

/* sort l[0]...l[n-1] */
void
sort (unsigned long *l, unsigned long n)
{
  unsigned long m, i, k;
  unsigned long *l1, *l2;

  if (n <= 1)
    return;

  m = n / 2;
  sort (l, m);
  sort (l2 = l + m, k = n - m);

  /* copy first sublist in l1[] */
  l1 = (unsigned long*) malloc (m * sizeof (unsigned long));
  for (i = 0; i < m; i++)
    l1[i] = l[i];

  /* now merge l1:m and l2:k */
  i = 0;
  while (m > 0 || k > 0)
    {
      if (m > 0 && k > 0)
	{
	  if (l1[0] <= l2[0])
	    {
	      l[i++] = l1[0];
	      l1 ++;
	      m --;
	    }
	  else
	    {
	      l[i++] = l2[0];
	      l2 ++;
	      k --;
	    }
	}
      else if (m > 0)
	{
	  l[i++] = l1[0];
	  l1 ++;
	  m --;
	}
      else /* k > 0 */
	{
	  l[i++] = l2[0];
	  l2 ++;
	  k --;
	}
    }
  free (l1 - n / 2);
}

/* remove duplicates */
unsigned long
uniq (unsigned long *l, unsigned long n)
{
  unsigned long i, j;
  
  j = 0;
  for (i = 0; i < n; i++)
    {
      while (i + 1 < n && l[i] == l[i + 1])
	i ++;
      l[j++] = l[i];
    }
  return j;
}

void
init_relation (relation r)
{
  r->q = 0;
  mpz_init_set_ui (r->x, 1);
  r->a = 0;
  r->p = NULL;
  r->e = NULL;
}

void
clear_relation (relation r)
{
  mpz_clear (r->x);
  if (r->p != NULL)
    free (r->p);
  if (r->e != NULL)
    free (r->e);
}

void
swap_relation (relation r, relation s)
{
  int *t;
  unsigned int *u;
  unsigned int v;

  mpz_swap (r->x, s->x);
  t = r->p; r->p = s->p; s->p = t;
  u = r->e; r->e = s->e; s->e = u;
  v = r->a; r->a = s->a; s->a = v;
  v = r->n; r->n = s->n; s->n = v;
}

/* return EOF if end-of-file or error, 0 otherwise */
int
read_partial (relation r, FILE *fp, int with_q)
{
  unsigned long q;
  unsigned int n = 0, e;
  int p;

  if (with_q)
    {
      if (fscanf (fp, "%lu", &q) == EOF)
	return EOF;
      r->q = q;
    }
  if (mpz_inp_str (r->x, fp, 10) == 0)
    return EOF;

  while (!feof (fp))
    {
      if (fscanf (fp, "%d", &p) == EOF)
        return EOF;
      if (p == 0) /* end of line */
        return 0;
      if (fscanf (fp, "%u", &e) == EOF)
        return EOF;
      realloc_relation (r, n + 1);
      r->p[n] = p;
      r->e[n] = e;
      r->n = ++ n;
    }
  return EOF;
}

/* return EOF if end-of-file or error, 0 otherwise */
int
read_partial2 (relation r, FILE *fp, unsigned long *pp)
{
  unsigned long q;
  unsigned int n = 0, e;
  int p;

  if (fscanf (fp, "%lu", pp) == EOF)
    return EOF;
  if (fscanf (fp, "%lu", &q) == EOF)
    return EOF;
  r->q = q;
  if (mpz_inp_str (r->x, fp, 10) == 0)
    return EOF;

  while (!feof (fp))
    {
      if (fscanf (fp, "%d", &p) == EOF)
        return EOF;
      if (p == 0) /* end of line */
        return 0;
      if (fscanf (fp, "%u", &e) == EOF)
        return EOF;
      realloc_relation (r, n + 1);
      r->p[n] = p;
      r->e[n] = e;
      r->n = ++ n;
    }
  return EOF;
}

/* return EOF if end-of-file or error, 0 otherwise */
int
read_full (relation r, FILE *fp, unsigned long fbn)
{
  unsigned int n = 0, e;
  int p;

  r->q = 0;
  if (mpz_inp_str (r->x, fp, 10) == 0)
    return EOF;

  while (!feof (fp))
    {
      if (fscanf (fp, "%d", &p) == EOF)
        return EOF;
      if (p == 0) /* end of line */
        return 0;
      /* check prime index */
      if (p != -1 && (p <= 0 || fbn < (unsigned long) p))
        {
          fprintf (stderr, "Error, wrong prime index in full: %d\n", p);
          exit (1);
        }
      if (fscanf (fp, "%u", &e) == EOF)
        return EOF;
      realloc_relation (r, n + 1);
      r->p[n] = p;
      r->e[n] = e;
      r->n = ++ n;
    }
  return EOF;
}

void
realloc_relation (relation r, unsigned int n)
{
  if (r->a < n)
    {
      r->a = n;
      r->p = (int*) realloc (r->p, n * sizeof (int));
      if (r->p == NULL)
        {
          fprintf (stderr, "Error, no more memory available\n");
          exit (1);
        }
      r->e = (unsigned int*) realloc (r->e, n * sizeof (unsigned int));
      if (r->e == NULL)
        {
          fprintf (stderr, "Error, no more memory available\n");
          exit (1);
        }
    }
}

void
copy_relation (relation r, relation s)
{
  unsigned int i;

  r->q = s->q;
  mpz_set (r->x, s->x);
  r->n = s->n;
  realloc_relation (r, s->n);
  for (i = 0; i < s->n; i++)
    {
      r->p[i] = s->p[i];
      r->e[i] = s->e[i];
    }
}

/* replace q, x by 1, x/q mod N */
void
normalize_relation (relation r, mpz_t N)
{
  mpz_t t;
  mpz_init (t);
  mpz_set_ui (t, r->q);
  mpz_invert (t, t, N);
  mpz_mul (t, t, r->x);
  mpz_mod (r->x, t, N);
  r->q = 1;
  mpz_clear (t);
}

void
output_relation (FILE *fp, relation r)
{
  unsigned int i;

  fprintf (fp, "%lu ", r->q);
  mpz_out_str (fp, 10, r->x);
  for (i = 0; i < r->n; i++)
    fprintf (fp, " %d %u", r->p[i], r->e[i]);
  fprintf (fp, " 0\n");
}

void
output_relation_without_q (FILE *fp, relation r)
{
  unsigned int i;

  mpz_out_str (fp, 10, r->x);
  for (i = 0; i < r->n; i++)
    fprintf (fp, " %d %u", r->p[i], r->e[i]);
  fprintf (fp, " 0\n");
}

void
output_relation2 (FILE *fp, relation r, mpz_t q, mpz_t N)
{
  unsigned int i;

  mpz_invert (q, q, N);
  mpz_mul (q, q, r->x);
  mpz_mod (r->x, q, N);
  mpz_set_ui (q, 1);
  mpz_out_str (fp, 10, r->x);
  for (i = 0; i < r->n; i++)
    fprintf (fp, " %d %u", r->p[i], r->e[i]);
  fprintf (fp, " 0\n");
}

/* for each prime index fb[0..fbn], count the number of occurrences of 
   an odd exponent in the relations */
void
count_odd_exponents (int *fb, unsigned long fbn, relation *rels,
		     unsigned long nrels)
{
  unsigned long i, j;
  
  for (i = 0; i <= fbn; i++)
    fb[i] = 0;
  for (i = 0; i < nrels; i++)
    for (j = 0; j < rels[i]->n; j++)
      if (rels[i]->e[j] % 2)
        fb[INDEX(rels[i]->p[j])] ++;
}

/* read number to factor from 'params' file */
void
read_params (mpz_t N, unsigned long *B, unsigned long *LP, char *params)
{
  FILE *fp;
  mpz_t lp, bb;

  fp = fopen (params, "r");
  if (mpz_inp_str (N, fp, 10) == 0)
    {
      fprintf (stderr, "Error, cannot read number to factor\n");
      exit (1);
    }
  if (fscanf (fp, "%lu", B) == EOF)
    {
      fprintf (stderr, "Error, cannot read factor base bound\n");
      exit (1);
    }
  mpz_init (lp);
  mpz_init (bb);
  if (mpz_inp_str (lp, fp, 10) == 0)
    {
      mpz_set_ui (lp, *B);
      mpz_mul_ui (lp, lp, 128); 
    }
  /* LP must fit in a "unsigned long" (so that large primes can be
     represented by this type), and be less than B^2 (so that when
     a residue is less than LP, it is necessarily prime */
  mpz_ui_pow_ui (bb, *B, 2);
  if (mpz_fits_ulong_p (lp) == 0 || mpz_cmp (lp, bb) > 0)
    {
      fprintf (stderr, "Error, large prime bound is too large\n");
      exit (1);
    }
  *LP = mpz_get_ui (lp);
  mpz_clear (lp);
  mpz_clear (bb);
  fclose (fp);
}

/* multiply s and t (mod N) into r, ignoring the large prime (if any).
   Overlap is not allowed between r and (s or t).
*/
void
multiply_relation (relation r, relation s, relation t, mpz_t N)
{
  unsigned int sn = s->n, tn = t->n, si = 0, ti = 0, ri = 0;

  if (r == s || r == t)
    abort ();

  mpz_mul (r->x, s->x, t->x);
  mpz_mod (r->x, r->x, N);
  realloc_relation (r, sn + tn);
  while (si < sn || ti < tn)
    {
      if (si < sn && ti < tn)
        {
          if (s->p[si] < t->p[ti])
            {
              r->p[ri] = s->p[si];
              r->e[ri] = s->e[si];
              si ++;
            }
          else if (t->p[ti] < s->p[si])
            {
              r->p[ri] = t->p[ti];
              r->e[ri] = t->e[ti];
              ti ++;
            }
          else /* s->p[si] = t->p[ti] */
            {
              r->p[ri] = s->p[si];
              r->e[ri] = s->e[si] + t->e[ti];
              si ++;
              ti ++;
            }
          ri ++;
        }
      else if (si < sn) /* t is empty */
        {
          r->p[ri] = s->p[si];
          r->e[ri] = s->e[si];
          si ++;
          ri ++;
        }
      else /* s empty */
        {
          r->p[ri] = t->p[ti];
          r->e[ri] = t->e[ti];
          ti ++;
          ri ++;
        }
    }
  r->n = ri;
}

/* put in fb[1..fbn] the factor base primes from file 'fbase',
   starting by index 1 (by convention, fb[0] represents -1).
   Return the number of primes in 'fbase' (-1 not counted).
*/
unsigned long
read_factor_base (int *fb, char *fbase, unsigned long fbn)
{
  FILE *fp;
  unsigned long n = 0;
  int p;

  fb[0] = -1; /* implicit */
  fp = fopen (fbase, "r");
  if (fp == NULL)
    {
      fprintf (stderr, "Error, cannot open %s\n", fbase);
      exit (1);
    }
  while (fscanf (fp, "%d", &p) != EOF)
    {
      n ++;
      if (n > fbn)
        {
          fprintf (stderr, "Error, too many primes in factor base\n");
          fprintf (stderr, "fbn=%lu n=%lu p=%d\n", fbn, n, p);
          exit (1);
        }
      fb[n] = p;
    }
  fclose (fp);
  return n;
}

/* check relation is correct, i.e. x^2 = p1^e1*p2^e2*... mod N.
   Return 0 iff the relation is ok.
   (q is not used here.)
*/
int
check_relation (relation r, unsigned int *fb, mpz_t N)
{
  mpz_t X, Y;
  unsigned long i, j;
  int p, res;

  mpz_init (X);
  mpz_init (Y);
  if (r->q == 0) /* full */
    mpz_set_ui (Y, 1);
  else /* partial */
    {
      mpz_set_ui (Y, r->q);
      mpz_mul_ui (Y, Y, r->q);
      mpz_mod (Y, Y, N);
    }
  mpz_mul (X, r->x, r->x);
  mpz_mod (X, X, N);
  for (i = 0; i < r->n; i++)
    for (j = 0; j < r->e[i]; j++)
      {
	p = (r->p[i] == -1) ? -1 : (int) fb[r->p[i]];
	mpz_mul_si (Y, Y, p);
	mpz_mod (Y, Y, N);
      }
  res = mpz_cmp (X, Y);
  mpz_clear (X);
  mpz_clear (Y);
  return res;
}

int
cputime (void)
{
  struct rusage rus;

  getrusage (0, &rus);
  return rus.ru_utime.tv_sec * 1000 + rus.ru_utime.tv_usec / 1000;
}

unsigned long
nextprime (unsigned long n)
{
  mpz_t p;

  mpz_init_set_ui (p, n);
  mpz_nextprime (p, p);
  n = mpz_get_ui (p);
  mpz_clear (p);
  return n;
}

trie_t*
insert_prime (trie_t *t, unsigned long p)
{
  if (t == NULL)
    {
      t = (trie_t*) malloc (sizeof (trie_t));
      t->p = p;
      t->left = NULL;
      t->right = NULL;
    }
  else if (p < t->p)
    t->left = insert_prime (t->left, p);
  else if (p > t->p)
    t->right = insert_prime (t->right, p);
  else /* p = t->p: do nothing */
    {
    }
  return t;
}

unsigned long
trie_size (trie_t *t)
{
  if (t == NULL)
    return 0;
  else
    return 1 + trie_size (t->left) + trie_size (t->right);
}

void
trie_clear (trie_t *t)
{
  if (t != NULL)
    {
      trie_clear (t->left);
      trie_clear (t->right);
      free (t);
    }
}

int
trie_contains (trie_t *t, unsigned long p)
{
  if (t == NULL)
    return 0;
  else if (p == t->p)
    return 1;
  else if (p < t->p)
    return trie_contains (t->left, p);
  else
    return trie_contains (t->right, p);
}
