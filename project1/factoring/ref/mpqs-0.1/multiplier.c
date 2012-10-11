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

/* evaluate the Knuth-Schroeppel function, cf Robert D. Silverman,
   "The Multiple Polynomial Quadratic Sieve", Math. of Comp. volume 48,
    number 177, 1987, page 335 */
unsigned long
find_multiplier (mpz_t N, double B)
{
  unsigned long k, bestk = 1;
  double p, f, g, maxf = 0.0;
  mpz_t kN;
  
  mpz_init (kN);
  for (k = 1; k < 100; k = nextprime (k))
    {
      mpz_mul_ui (kN, N, k);
      /* FIXME: Silverman writes "if N = 1 mod 8" but isn't it kN instead? */
      if (mpz_kronecker_ui (kN, 2) == 1 && mpz_fdiv_ui (kN, 8) == 1)
        f = 2.0 * log (2.0);
      else
        f = 0.0;
      for (p = getprime (2.0); p <= B; p = getprime (p))
        {
          if (mpz_kronecker_ui (kN, (unsigned long) p) == 1)
            {
              g = ((k % (unsigned long) p) == 0) ? (1.0 / p) : (2.0 / p);
              f += g * log (p);
            }
        }
      f -= 0.5 * log ((double) k);
      if (f > maxf)
        {
          maxf = f;
          bestk = k;
        }
      getprime (0.0); /* free prime buffer */
    }
  mpz_clear (kN);
  
  return bestk;
}
