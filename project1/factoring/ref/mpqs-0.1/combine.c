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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h> /* for ULONG_MAX */
#include <math.h>
#include "mpqs.h"

#define MAXSPECIAL 2147483648UL /* maximal number of special nodes */
#define IS_SPECIAL(x) (x > (ULONG_MAX - MAXSPECIAL))

unsigned long special = 0; /* number of special nodes */
unsigned long components = 0; /* number of connected components <= special */

static void
usage (void)
{
  fprintf (stderr, "Usage: combine2 -p <params>\n");
  exit (1);
}

/* read one partial relation, and put the large prime in p[0].
   Return EOF if end of file.
 */
static int
readline1 (FILE *fp, unsigned long *p)
{
  if (fscanf (fp, "%lu", p) == EOF)
    return EOF;
  while (fgetc (fp) != '\n');
  return 1;
}

/* read one pp, and put the large primes in p[0] and p[1] */
static int
readline2 (FILE *fp, unsigned long *p, unsigned long *q)
{
  if (fscanf (fp, "%lu %lu", p, q) == EOF)
    return EOF;
  while (fgetc (fp) != '\n');
  return 1;
}

static trie_t*
select_primes (char *par, char *pps, unsigned long LP)
{
  unsigned long n;
  FILE *fp;
  unsigned long p, q, h;
  mpz_t T;
  unsigned long nfp = 0, npp = 0;
  trie_t *t = NULL;
  
  /* choose hash size */
  n = (unsigned long) ((double) LP / log ((double) LP));
  n = nextprime (n);

  mpz_init (T);

  fp = fopen (par, "r");
  while (!feof (fp))
    {
      if (readline1 (fp, &p) == EOF)
	break;
      nfp ++;
      h = p % n;
      if (mpz_tstbit (T, h))
	t = insert_prime (t, p);
      else
	mpz_setbit (T, h);
    }
  fclose (fp);

  fp = fopen (pps, "r");
  while (!feof (fp))
    {
      if (readline2 (fp, &p, &q) == EOF)
	break;
      npp ++;
      h = p % n;
      if (mpz_tstbit (T, h))
	t = insert_prime (t, p);
      else
	mpz_setbit (T, h);
      h = q % n;
      if (mpz_tstbit (T, h))
	t = insert_prime (t, q);
      else
	mpz_setbit (T, h);
    }
  fclose (fp);

  fprintf (stderr, "read %lu fps, %lu pps\n", nfp, npp);
  fflush (stderr);

  mpz_clear (T);

  return t;
}

/* return index i of p in H: either p is already in H, then H[i]=p,
   or p is not, and H[i]=0 */
static unsigned long
seek (relation *H, unsigned long n, unsigned long p)
{
  unsigned long h;

  h = p % n;
  while (H[h]->q != p && H[h]->q != 0)
    h = (h + 1) % n;
  return h;
}

/* root nodes are those such that F[i] is special */
static unsigned long
root (unsigned long *F, unsigned long i)
{
  while (IS_SPECIAL(F[i]) == 0)
    i = F[i];
  return i;
}

#ifdef DEBUG
static unsigned long
print_fathers (unsigned long *F, relation *H, unsigned long i)
{
  while (IS_SPECIAL(F[i]) == 0)
    {
      printf ("%lu->", H[i]->q);
      i = F[i];
    }
  printf ("%lu\n", H[i]->q);
  return i;
}
#endif

/* return 1 if a new cycle is found, 0 otherwise */
static int
insert (relation *H, unsigned long *F, unsigned long n,
        unsigned long p, unsigned long q, relation r, mpz_t N, FILE *out)
{
  unsigned long hp, hq;

  hp = seek (H, n, p);
  hq = seek (H, n, q);
  if (hp == hq) /* necessarily H[hp]->q = H[hq]->q = 0 */
    {
      do
        {
          hq ++;
        }
      while (H[hq]->q != 0);
    }
  if (H[hp]->q == p && H[hq]->q == q) /* both p and q are already in */
    {
      if (root (F, hp) == root (F, hq))
	{
	  relation s;
	  mpz_t Q;
	  init_relation (s);
	  mpz_init_set_ui (Q, 1); /* accumulates large primes */
	  /* first accumulate in r all relations from hp to root */
	  while (IS_SPECIAL(F[hp]) == 0)
	    {
	      copy_relation (s, r);
	      multiply_relation (r, s, H[hp], N);
	      mpz_mul_ui (Q, Q, H[hp]->q);
	      hp = F[hp];
	    }
	  /* special nodes have no relation, only a large prime */
	  mpz_mul_ui (Q, Q, H[hp]->q);
	  /* then from hq to root */
	  while (IS_SPECIAL(F[hq]) == 0)
	    {
	      copy_relation (s, r);
	      multiply_relation (r, s, H[hq], N);
	      mpz_mul_ui (Q, Q, H[hq]->q);
	      hq = F[hq];
	    }
          /* H[hp]->q = H[hq]->q, we want to accumulate it only once */
	  output_relation2 (out, r, Q, N);
	  clear_relation (s);
	  mpz_clear (Q);
	  return 1;
	}
      else /* join the component of q to that of p */
        {
          unsigned long t, hr = F[hq]; /* father of q */

	  components --;
          F[hq] = hp;
          while (IS_SPECIAL(hr) == 0)
            {
              swap_relation (r, H[hq]);
              t = F[hr];
              F[hr] = hq;
              hq = hr;
              hr = t;
            }
          swap_relation (r, H[hq]);
        }
    }
  else if (H[hp]->q == p) /* q is new: join p-q to component of p */
    {
      copy_relation (H[hq], r);
      H[hq]->q = q;
      F[hq] = hp;
    }
  else if (H[hq]->q == q) /* p is new: join p-q to component of q */
    {
      copy_relation (H[hp], r);
      H[hp]->q = p;
      F[hp] = hq;
    }
  else /* both p and q are new: create a new special node */
    {
      components ++;
      H[hp]->q = p;
      copy_relation (H[hq], r);
      H[hq]->q = q;
      F[hp] = ULONG_MAX - special; /* hp is special */
      F[hq] = hp;
      special ++;
    }
  return 0;
}

static void
second_pass (char *par, char *pps, trie_t *t, mpz_t N, char *out)
{
  FILE *fp, *fp2;
  unsigned long nfp = 0, npp = 0;
  unsigned long p, q;
  relation *H;
  unsigned long *F;
  unsigned long n, i;
  unsigned long cycles = 0, cycles_fp;
  relation r;

  init_relation (r);

  n = trie_size (t);
  fprintf (stderr, "found %lu large primes appearing at least twice\n", n);

  n = nextprime (n + n / 2);

  H = (relation*) malloc (n * sizeof (relation));
  F = (unsigned long*) malloc (n * sizeof (unsigned long));

  /* init table */
  for (i = 0; i < n; i++)
    {
      init_relation (H[i]);
      F[i] = 0;
    }

  fp = fopen (par, "r");
  fp2 = fopen (out, "w");
  while (!feof (fp))
    {
      if (read_partial (r, fp, 1) == EOF)
	break;
      if (trie_contains (t, r->q))
	{
	  nfp ++;
	  /* for full-partials, the 1st prime is 1 */
	  cycles += insert (H, F, n, 1, r->q, r, N, fp2);
	}
    }
  fclose (fp);
  fprintf (stderr, "found %lu cycles from fps\n", cycles_fp = cycles);
  fflush (stderr);

  fp = fopen (pps, "r");
  while (!feof (fp))
    {
      if (read_partial2 (r, fp, &p) == EOF)
	break;
      q = r->q;
      if (trie_contains (t, p) && trie_contains (t, q))
	{
	  npp ++;
	  cycles += insert (H, F, n, p, q, r, N, fp2);
	}
    }
  fclose (fp);

  fclose (fp2);
  fprintf (stderr, "kept %lu fps and %lu pps\n", nfp, npp);
  fprintf (stderr, "%lu connected component(s)\n", components);
  fprintf (stderr, "found %lu cycles from pps, total %lu\n",
	   cycles - cycles_fp, cycles);
  fflush (stderr);

  /* clear table */
  for (i = 0; i < n; i++)
    clear_relation (H[i]);
  clear_relation (r);

  free (H);
  free (F);
}

int
main (int argc, char *argv[])
{
  char *par = NULL;
  char *pps = NULL;
  char *output = NULL;
  char *params = NULL;
  int st = cputime ();
  unsigned long B;  /* factor base bound */
  unsigned long LP; /* large prime bound */
  mpz_t N;          /* number to factor */
  trie_t *t;        /* trie of primes appearing at least twice in fps or pps */

  while (argc > 1)
    {
      if (argc > 2 && strcmp (argv[1], "-p") == 0)
        {
          params = argv[2];
          argv += 2;
          argc -= 2;
        }
      else
        usage ();
    }

  if (params == NULL)
    usage ();

  mpz_init (N);

  par = (char*) malloc (strlen (params) + strlen (PARTIALS) + 1);
  strcpy (par, params);
  strcat (par, PARTIALS);
  pps = (char*) malloc (strlen (params) + strlen (PPS) + 1);
  strcpy (pps, params);
  strcat (pps, PPS);
  output = (char*) malloc (strlen (params) + strlen (CYCLES) + 1);
  strcpy (output, params);
  strcat (output, CYCLES);

  read_params (N, &B, &LP, params);

  t = select_primes (par, pps, LP);

  second_pass (par, pps, t, N, output);

  fprintf (stderr, "combining partials took %1.1fs\n", 
	   (double) (cputime () - st) / 1e3);

  free (pps);
  free (output);

  mpz_clear (N);

  return 0;
}
