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

static void
usage (void)
{
  fprintf (stderr, "Usage: buildmatrix -p <params> [-keep <nnn>]\n");
  exit (1);
}

static unsigned long
weight (relation r, unsigned int *fb)
{
  unsigned long w = 0, j;

  for (j = 0; j < r->n; j++)
    w += (r->e[j] % 2) && (fb[INDEX(r->p[j])] >= 1);
  return w;
}

int
main (int argc, char *argv[])
{
  char *params = NULL, *fbase, *fulls, *cycles, *matrix;
  unsigned long fbn; /* number of primes in factor base */
  unsigned long ncycles; /* number of cycles */
  FILE *fp;
  relation *rels;
  unsigned long nfulls;
  unsigned long nrels, nrels0;
  unsigned long i, j, k;
  unsigned int *fb;
  unsigned long zero, discard, nprimes;
  long keep = 10, excess;
  mpz_t N;
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
      else if (argc > 2 && strcmp (argv[1], "-keep") == 0)
        {
          keep = atoi (argv[2]);
          argv += 2;
          argc -= 2;
        }
      else
        usage ();
    }

  if (params == NULL)
    usage ();

  mpz_init (N);

  read_params (N, &B, &LP, params);

  fbase = (char*) malloc (strlen (params) + strlen (FBASE) + 1);
  strcpy (fbase, params);
  strcat (fbase, FBASE);

  fulls = (char*) malloc (strlen (params) + strlen (FULLS) + 1);
  strcpy (fulls, params);
  strcat (fulls, FULLS);

  cycles = (char*) malloc (strlen (params) + strlen (CYCLES) + 1);
  strcpy (cycles, params);
  strcat (cycles, CYCLES);

  matrix = (char*) malloc (strlen (params) + strlen (MATRIX) + 1);
  strcpy (matrix, params);
  strcat (matrix, MATRIX);

  fbn = count_lines (fbase);
  fprintf (stderr, "factor base has %lu primes\n", fbn);
  fb = (unsigned int*) malloc ((fbn + 1) * sizeof (unsigned int));
  read_factor_base ((int*) fb, fbase, fbn);

  ncycles = count_lines (cycles);
  fprintf (stderr, "%lu cycles\n", ncycles);

  nfulls = count_lines (fulls);
  fprintf (stderr, "%lu full relations\n", nfulls);

  nrels0 = nrels = nfulls + ncycles;
  fprintf (stderr, "%lu total relations\n", nrels0);

  rels = (relation*) malloc (nrels0 * sizeof (relation));
  for (i = 0; i < nrels0; i++)
    init_relation (rels[i]);

  /* first nfulls relations are fulls, next ncycles are cycles */
  /* first fbn primes are from factor base */

  /* read fulls */
  fp = fopen (fulls, "r");
  for (i = 0; i < nfulls; i++)
    {
      read_full (rels[i], fp, fbn);
      if (check_relation (rels[i], fb, N))
	{
	  fprintf (stderr, "Error, wrong full relation:\n");
          output_relation_without_q (stderr, rels[i]);
	  exit (1);
	}
    }
  fclose (fp);

  /* read cycles */
  fp = fopen (cycles, "r");
  for (i = 0; i < ncycles; i++)
    {
      j = nfulls + i;
      read_partial (rels[j], fp, 0);
      if (check_relation (rels[j], fb, N))
        {
          fprintf (stderr, "Error, wrong cycle relation:\n");
          output_relation_without_q (stderr, rels[j]);
          /* exit (1); */
          i --;
          ncycles --;
        }
    }
  fclose (fp);

  /* loop to count number of odd exponents for each factor base prime */
  nprimes = fbn + 1;
  zero = 0;
  excess = (long) nrels - (long) nprimes;
  do
    {
      fprintf (stderr, "remains %lu rels on %lu primes (excess %ld)\n", nrels,
	       nprimes - zero, excess);
      /* FIXME: avoid recomputing each time, by updating wrt removed
	 relations */
      count_odd_exponents ((int*) fb, fbn, rels, nrels);

      /* count "zero" primes */
      zero = 0;
      for (i = 0; i <= fbn; i++)
	zero += fb[i] <= 1;

      /* remove relations with "single" primes */
      k = 0;
      for (i = 0; i < nrels; i++)
	{
	  for (j = 0; j < rels[i]->n; j++)
	    if (rels[i]->e[j] % 2 && fb[INDEX(rels[i]->p[j])] == 1)
              break;
	  if (j == rels[i]->n)
	    {
	      copy_relation (rels[k], rels[i]);
	      k ++;
	    }
	}
      discard = nrels - k;
      if (discard > 0)
	fprintf (stderr, "discarded %lu singleton(s)\n", discard);
      nrels = k;

      excess = (long) nrels - (long) (nprimes - zero);
      if (discard == 0 && excess > keep)
	{ /* remove (excess-keep) heavier relations */
	  unsigned long w, minw, nmin;
	  unsigned long *remove, nremove = 0, toremove = excess - keep;
	  remove = (unsigned long*) malloc (toremove * sizeof (unsigned long));
	  minw = ULONG_MAX;
	  nmin = 0;
	  for (i = 0; i < nrels; i++)
	    {
	      w = weight (rels[i], fb);
	      if (nremove < toremove)
		{
		  if (w < minw)
		    {
		      minw = w; /* min, weight of current removed relations */
		      nmin = 1; /* number of min. weight relations */
		    }
		  else if (w == minw)
		    nmin ++;
		  remove[nremove++] = i;
		}
	      else if (w > minw) /* compare with current heavy relations */
		{
		  /* find a relation of minimum weight */
		  for (j = 0; weight (rels[remove[j]], fb) != minw; j++);
		  /* replace it by the current one */
		  remove[j] = i;
		  if (-- nmin == 0) /* update minw, nmin */
		    {
		      minw = w + 1; /* surely larger than the minimum */
		      for (j = 0; j < nremove; j++)
			if ((w = weight (rels[remove[j]], fb)) < minw)
			  {
			    minw = w;
			    nmin = 1;
			  }
			else if (w == minw)
			  nmin ++;
		    }
		}
	    }
	  sort (remove, nremove);
	  for (j = 1; j <= nremove; j--)
	    { /* swap rels[remove[nremove-j]] and rels[nrels-j] */
	      if (remove[nremove - j] < nrels - j)
		copy_relation (rels[remove[nremove - j]], rels[nrels - j]);
	    }
	  free (remove);
	  fprintf (stderr, "removed %lu heavy relations\n", nremove);
	  discard = nremove;
	  nrels -= nremove;
	  excess -= nremove;
	}
    }
  while (discard > 0 || excess > keep);
  free (fb);

  fp = fopen (matrix, "w");
  /* output relations */
  for (i = 0; i < nrels; i++)
    output_relation_without_q (fp, rels[i]);
  fclose (fp);

  for (i = 0; i < nrels0; i++)
    clear_relation (rels[i]);
  free (rels);

  free (fbase);
  free (fulls);
  free (cycles);
  free (matrix);

  mpz_clear (N);

  fprintf (stderr, "building matrix took %1.1fs\n",
           (double) (cputime () - st) / 1e3);

  return 0;
}
