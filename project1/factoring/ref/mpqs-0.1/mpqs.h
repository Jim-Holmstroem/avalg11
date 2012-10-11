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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* for memset */
#include <limits.h> /* ULONG_MAX */
#include <sys/types.h>
#include <sys/resource.h>
#include "gmp.h"

#define FBASE    ".fbase"
#define FULLS    ".fulls"
#define PARTIALS ".partials"
#define PPS      ".pps"
#define CYCLES   ".cycles"
#define AVAL     ".a_val"
#define STATUS   ".status"
#define ADATA    ".adata"
#define MATRIX   ".matrix"

#define INDEX(p) ((p == -1) ? 0 : p)

typedef struct {
  unsigned long q; /* large prime */
  mpz_t x;
  unsigned int n; /* number of primes stored */
  unsigned int a; /* number of primes allocated */
  int *p; /* primes and/or indices */
  unsigned int *e;  /* exponents */
} relation_struct;
typedef relation_struct relation[1];

typedef struct _trie_t {
  unsigned long p; /* large prime */
  struct _trie_t *left;
  struct _trie_t *right;
} trie_t;

unsigned long count_lines (char*);
void sort (unsigned long*, unsigned long);
unsigned long uniq (unsigned long*, unsigned long);
void init_relation (relation);
void clear_relation (relation);
int read_partial (relation, FILE*, int);
int read_partial2 (relation, FILE *, unsigned long *);
int read_full (relation, FILE*, unsigned long);
void realloc_relation (relation, unsigned int);
void copy_relation (relation, relation);
void output_relation (FILE *, relation);
void output_relation_without_q (FILE *, relation);
void output_relation2 (FILE *, relation, mpz_t, mpz_t);
void count_odd_exponents (int*, unsigned long, relation*, unsigned long);
void read_params (mpz_t, unsigned long*, unsigned long*, char*);
void multiply_relation (relation, relation, relation, mpz_t);
unsigned long read_factor_base (int *, char *, unsigned long);
int check_relation (relation, unsigned int *, mpz_t);
int cputime (void);
unsigned long nextprime (unsigned long);
trie_t* insert_prime (trie_t *, unsigned long);
unsigned long trie_size (trie_t *);
void trie_clear (trie_t *);
int trie_contains (trie_t *, unsigned long);
void normalize_relation (relation, mpz_t);
void swap_relation (relation, relation);
double getprime (double);
unsigned long find_multiplier (mpz_t, double);
