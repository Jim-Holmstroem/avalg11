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
#include "mpqs.h"

static void
usage (void)
{
  fprintf (stderr, "Usage: gauss -p <params>\n");
  exit (1);
}

/* Perform naive Gaussian elimination on matrix M (i rows, j columns).
   Each row is represented by a bit-vector: M[i,j] is bit (j-1) of M[i].
   Return the row of the first dependency */
static int
gauss (mpz_t *M, unsigned long i, unsigned long j)
{
  unsigned long k, l, m;
  unsigned long p = 1;

#define WIDTH 80

  k = l = 0;
  /* current row is k, column is l */
  while (k < i && l < j)
    {
      if (k == (p * i) / WIDTH)
	{
	  fprintf (stderr, ".");
	  fflush (stderr);
	  p ++;
	}
      /* search pivot in column l */
      for (m = k; m < i && mpz_tstbit (M[m], l) == 0; m++);
      if (m < i) /* found pivot */
        {
          mpz_swap (M[k], M[m]);
          for (m = k + 1; m < i; m++)
            if (mpz_tstbit (M[m], l))
              mpz_xor (M[m], M[m], M[k]);
          k++;
        }
      /* if no pivot found, goto next column */
      l++;
    }
  fprintf (stderr, "\n");
  fflush (stderr);
  return k;
}

/* return non-zero iff a non-trivial factor was found */
static int
treat_dependency (mpz_t M, mpz_t N, relation *rels, unsigned long ncols,
		  int *fb)
{
  relation r, s;
  unsigned long i, j;
  int first = 1;
  mpz_t Q; /* collects the large primes */
  mpz_t g;
  int p, res;

  init_relation (r);
  init_relation (s);
  mpz_init_set_ui (Q, 1);
  for (j = ncols; j <= mpz_sizeinbase (M, 2); j++)
    if (mpz_tstbit (M, j)) /* multiply by relation (j-ncols) */
      {
	i = j - ncols;
	if (first)
	  {
	    copy_relation (r, rels[i]);
	    first = 0;
	  }
	else
	  {
	    copy_relation (s, r);
	    multiply_relation (r, s, rels[i], N);
	  }
	if (rels[i]->q != 0)
	  {
	    mpz_mul_ui (Q, Q, rels[i]->q);
	    mpz_mod (Q, Q, N);
	  }
      }

  /* check exponents are all even, and take square root */
  for (i = r->n; i-- > 0;)
    {
      j = r->e[i];
      p = r->p[i];
      if (j % 2)
	{
	  fprintf (stderr, "Error: odd exponent %lu for prime %d\n", j, p);
	  exit (1);
	  goto end;
	}
      j /= 2;
      while (j--)
	{
	  mpz_mul_si (Q, Q, fb[INDEX(p)]);
	  mpz_mod (Q, Q, N);
	}
    }

  mpz_init (g);
  mpz_add (g, r->x, Q);
  mpz_gcd (g, g, N);
  res = mpz_cmp_ui (g, 1) && mpz_cmp (g, N);
  fprintf (stderr, "gcd(X+Y,N)=");
  mpz_out_str (stderr, 10, g);
  fprintf (stderr, "\n");
  mpz_sub (g, r->x, Q);
  mpz_gcd (g, g, N);
  res = res || (mpz_cmp_ui (g, 1) && mpz_cmp (g, N));
  fprintf (stderr, "gcd(X-Y,N)=");
  mpz_out_str (stderr, 10, g);
  fprintf (stderr, "\n");
  mpz_clear (g);

 end:
  clear_relation (r);
  clear_relation (s);
  mpz_clear (Q);
  
  return res;
}

int
main (int argc, char *argv[])
{
  char *params = NULL, *matrix, *fbase;
  unsigned long nrows, ncols, i, j;
  relation *rels;
  FILE *fp;
  int *fb;
  unsigned int *fbi;
  unsigned long fbn = 0;
  mpz_t *M; /* bit matrix */
  mpz_t N; /* number to factor */
  int youpi = 0;
  unsigned long B, LP;
  int st = cputime ();

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

  fbase = (char*) malloc (strlen (params) + strlen (FBASE) + 1);
  strcpy (fbase, params);
  strcat (fbase, FBASE);

  matrix = (char*) malloc (strlen (params) + strlen (MATRIX) + 1);
  strcpy (matrix, params);
  strcat (matrix, MATRIX);

  nrows = count_lines (matrix);
  fprintf (stderr, "%lu relations\n", nrows);

  if (nrows == 0)
    return 0;

  fbn = count_lines (fbase);
  fb = (int*) malloc ((fbn + 1) * sizeof (int));
  read_factor_base (fb, fbase, fbn);
 
  mpz_init (N);
  read_params (N, &B, &LP, params);

  rels = (relation*) malloc (nrows * sizeof (relation));
  fp = fopen (matrix, "r");
  for (i = 0; i < nrows; i++)
    {
      init_relation (rels[i]);
      read_partial (rels[i], fp, 0);
      if (check_relation (rels[i], (unsigned int*) fb, N) != 0)
        {
          fprintf (stderr, "Error, wrong relation:\n");
          output_relation_without_q (stderr, rels[i]);
          exit (1);
        }
    }
  fclose (fp);

  /* find factor base primes with at least one odd exponent */
  count_odd_exponents (fb, fbn, rels, nrows);
  fbi = (unsigned int*) malloc ((fbn + 1) * sizeof (unsigned int));
  /* fbi[j] <= j is the column index of the j-th prime with odd exponent */
  for (i = 0, ncols = 0; i <= fbn; i++)
    {
      if (fb[i] == 1) /* we should have no singleton any more here */
	abort ();
      if (fb[i] >= 2)
	fbi[i] = ncols++;
    }
  fprintf (stderr, "%lu primes with odd exponent\n", ncols);

  M = (mpz_t*) malloc (nrows * sizeof (mpz_t));

  for (i = 0; i < nrows; i++)
    {
      int p;
      mpz_init (M[i]); /* sets to 0 */
      for (j = 0; j < rels[i]->n; j++)
	{
	  p = rels[i]->p[j];
	  if (fb[INDEX(p)] >= 2 && rels[i]->e[j] % 2) /* odd exponent */
	    mpz_setbit (M[i], fbi[INDEX(p)]);
	}
      /* add a bit to identify relation i */
      mpz_setbit (M[i], ncols + i);
    }

  i = gauss (M, nrows, ncols);
  /* i is the index of the first dependency */
  
  fprintf (stderr, "Found %lu dependencies\n", nrows - i);

  read_factor_base (fb, fbase, fbn);

  for (j = i; j < nrows && youpi == 0; j++)
    {
      fprintf (stderr, "Process dependency %lu:\n", j - i + 1);
      youpi = treat_dependency (M[j], N, rels, ncols, fb);
    }

  mpz_clear (N);
  for (i = 0; i < nrows; i++)
    mpz_clear (M[i]);
  free (M);
  free (fbi);
  free (fb);
  for (i = 0; i < nrows; i++)
    clear_relation (rels[i]);
  free (rels);

  free (fbase);
  free (matrix);

  fprintf (stderr, "Gaussian elimination took %1.1fs\n",
           (double) (cputime () - st) / 1e3);

  return 0;
}
