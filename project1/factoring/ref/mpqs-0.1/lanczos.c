#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpqs.h"


/* Change randomness (to get different solutions) by putting an integer
   in the seed file
 */
#define SEED_FILE    "SEED"

#define SIZEOFLONG    ((signed) sizeof(long))
#define WORD_SIZE     (8*SIZEOFLONG)
#define WORD_SIZEM    (WORD_SIZE-1)
#define RBITS         16    /* For randomness selection of bits */
#define BAAT          4     /* BAAT bits at a time, must be a power of 2 
                             * Used to speed up multiplying transpose of 
                             * block vector with itself (I think this is
                             * a Coppersmith trick), and product of a block
                             * vector with block.
                             *
                             * Recommendation: Take BAAT as 4 or 8 .
                             */



/* for output progress reports: */
#define WIDTH       20




#define PrintV(Vec) { \
    register int i,j; \
\
    for (i=0; i<n1; i++) { \
        for (j=WORD_SIZEM; j>=0; j--) \
            if ( (Vec[i]>>j)&1 ) \
                printf("1"); \
            else \
                printf("0"); \
        printf("\n"); \
    } \
}


#define PrintNN(Vec) { \
    register int i,j; \
\
    for (i=0; i< WORD_SIZE; i++) { \
        for (j=WORD_SIZEM; j>=0; j--) \
            if ( (Vec[i]>> j) & 1) \
                printf("1"); \
            else \
                printf("0"); \
        printf("\n"); \
    } \
}


#define nswapbits(xx,b1,b2) { \
    if (((xx>>b1)&1) != ((xx>>b2)&1)) \
        xx ^= ((1UL<<b1)|(1UL<<b2)); \
}

#define nswapbits2(xx,b1,yy,b2) { \
    if (((xx>>b1)&1) != ((yy>>b2)&1)) { \
        xx ^= (1UL<<b1); \
        yy ^= (1UL<<b2); \
    }\
}


#define nxorbits(xx, bi1, bi2) {  /* bit bi1 ^= bit bi2 */ \
    if ((xx>>bi2)&1 ) \
        xx ^= (1UL<<  bi1); \
}

#define nxorbits2(xx, bi1, yy, bi2) {  /* bit bi1 ^= bit bi2 */ \
    if ((yy>>bi2) & 1 ) \
        xx ^= (1UL<<  bi1); \
}


static void
GetY(Y, n2)
unsigned long *Y;
int n2;
{

    register int i,k;
    register unsigned long x;

    for (i=n2; i; ) {
        x=0;
        for (k=0; k<SIZEOFLONG/2; k++) {
            /* RBITS bits at a time */
            x<<=RBITS;
            x |= ((i+k+random())& ((1UL<<RBITS)-1));
        }
        Y[--i]=x;
    }
}



/* Search through  array  between  &array[low]  and  &array[high]  for
 * the element  target  ( array  is sorted).
 */
static int
binarysearch(array, low, high, target)
int *array;
int low, high, target;
{

    register int l = low, h = high, m = (l+h)>>1;

    while (l < h) {
        if (array[m] == target)
            return m;
        if (array[m] < target)
            l = m + 1;
        else h = m;
        m = (l + h)>>1;
    }
    if (array[m] > target && m > l)
        return m-1;
    return m;
}



static int**
ReadB( size_fb, rel_num, matrix_fn, n1, n2 )
int size_fb;
int **rel_num;
char *matrix_fn;
int *n1, *n2;
{
    register int i,j,k;
    int ind, d, a, kill_count, e;
    int *mem;
    int *Bsize;
    int *Bcount;
    int **B;
    FILE *fp;
    char buffer[4096];
    int *kill;
    int eq;
    int c;


    fp = fopen(matrix_fn,"r");
    if (fp==NULL) {
        printf("file '%s' does not exist\n",matrix_fn);
        exit(1);
    }



    printf("Reading B\n"); fflush(stdout);

    Bsize = (int *) malloc((size_fb) * sizeof(int));
    Bcount = (int *) malloc((size_fb) * sizeof(int));
    B = (int **) malloc((size_fb) * sizeof(int*));

    for (i=0; i < size_fb; i++) {
        Bsize[i] = 2;
        B[i] = (int *) malloc(Bsize[i] * sizeof(int) );
        /* Bsize  keeps track of the amount of memory allocated for
           row  B[i]  of the matrix.  Bcount[i]  keeps track of the number
           of elements currently being used in that row.
         */
        Bcount[i] = 0;
    }

    /*
     * primes correspond to rows, equations correspond to cols
     *
     * n1  will keep track of the largest index within the factor base
     * n2  will keep track of the relation number
     */
    for (*n1=0,*n2=0; ; (*n2)++) {

        /* skip past first elements (x) of data file */
        while ( (c=getc(fp)) != ' ' )
            if (c == EOF)
                goto done;

        /* now read in data from the relation */
        while(1) {
            /* ind  represents the index of a factor base element */
	  fscanf(fp,"%d",&ind);
	  if (ind == 0)
	    /* end of relation */
	    break;
	  if (ind == -1)
	    ind = 0;
	  fscanf(fp,"%d",&e);
	  if (e & 1) {
                if (ind >= *n1) {
                    *n1 = ind + 1;
                    if (*n1 > size_fb) {
                        printf("ERROR! n1=%d, size_fb=%d\n",*n1,size_fb);
                        exit(1);
                    }
                }

                if (Bcount[ind] == Bsize[ind]) {
                    /* We've reached the memory limit for row  ind .
                       Allocate extra memory and copy things over
                     */
		    Bsize[ind] <<= 1;
		    B[ind] = (int*) realloc (B[ind], Bsize[ind] * sizeof(int));
                }
                B[ind][Bcount[ind]] = *n2;
		Bcount[ind] ++;

            }  /* if (e&1) */

        } /* while( 1 ) */
        fgets(buffer,sizeof(buffer),fp);
    }
    printf ("1\n");

done:
    free(Bsize);

    (*rel_num) = (int *) malloc(*n2 * sizeof(int));
    for (i=0; i < *n2; ++i)
        (*rel_num)[i] = i;
    kill = (int *) malloc((size_fb) * sizeof(int));

    printf("Read in %d by %d matrix\n",*n1,*n2); fflush(stdout);

    /*
     * The matrix allocated more memory than it actually needs.
     * Free the extra memory.  Also remove singletons which are
     * relations that contain primes that occur nowhere else in
     * the matrix.  Actually, Paul Zimmermann's code should already
     * remove singletons, so don't worry about repeatedly applying
     * this procedure.
     */
    j = 0;
    kill_count = 0;
    for (i=0; i< *n1; i++) {
        d = Bcount[i] + 1;
        if (d > 2 ) {
            mem = ( int *) malloc(sizeof(int) *d);
            /* -1 marks the end of the relation.  */
            mem[--d] = -1;
            for (k=d-1; k>= 0; k--)
                mem[k] = B[i][k];
            free(B[i]);
            Bcount[j] = Bcount[i];
            B[j++] = mem;
        }
        else {
            /* Either the row is a singleton, or else the row is empty.
               We can delete these rows completely.  The  kill  array
               holds the relation numbers that will be deleted.
             */
            if (d==2) {
                kill[kill_count++] = B[i][0];
                (*rel_num)[B[i][0]]= -1;
            }
            free(B[i]);
        }
    }
    /* Make  *n1  hold the number of primes that occur at least once */
    *n1 = j;

    /* sort  kill  array.  (insertion sort) */
    for (i=0; i < kill_count-1; i++) {
        j = i+1;
        while (j && ((k=kill[j]) < kill[j-1])) {
            kill[j] = kill[j-1];
            kill[--j] = k;
        }
    }
    kill[kill_count] = -1;


    if (kill_count) {
        /* remove the singletons:
         * scan each row of  B  to see if it has one of the elements from
         * the  kill  array.  If so, then the element is removed.
         */
        for (i=0; i < *n1; i++) {
            k=0;
            for (j=0; ;) {
                ind = B[i][j];
                if (ind== -1)
                    break;
                k = binarysearch(kill, k, kill_count-1, ind);
                if (ind == kill[k]) {
                    for (a=j; ; a++) {
                        eq = B[i][a+1];
                        B[i][a] = eq;
                        if (eq== -1) break;
                    }
                    (Bcount[i])--;
                }
                else j++;
            }
        }
    }


    free(Bcount);
    free(kill);

    /* Recompute the number of relations involved (after singletons have
     * been removed).  Let  n2  be this number.  Use the  rel_num  array
     * to keep track of the original relation numbers.
     */
    j=0;
    for (i=0; i < *n2; i++) {
        if ( (*rel_num)[i] == -1)
            continue;
        else {
            (*rel_num)[j++] = i;
        }
    }
    *n2 = j;

    /* Renumber the elements of the matrix to minimize memory requirements.
     * More precisely, number the columns from 0 to  n2-1  where  n2  was
     * computed above.
     */
    e=0;
    for (i=0; i < *n1; i++) {
        k=0;
        for (j=0; ; j++) {
            ind = B[i][j];
            if (ind== -1)
                break;
            a = ind;
            if (a >= *n2)
                a= *n2-1;
            k= binarysearch(*rel_num,k,a,ind);
            B[i][j] = k;
            ++e;
        }
    }

    printf("%d by %d matrix\n",*n1,*n2); fflush(stdout);
    printf("%d nonzero entries\n",e);
    printf("avg of %d nonzero entries per column\n",(int) (e/(double)*n2 + .5));

    return B;
}



static void
AtimesV(B, in_V, out_V, n1, n2)
int *B[];
unsigned long *in_V;
unsigned long *out_V;
int n1, n2;
{

    register int i,l;
    register unsigned long x;
    register int *ptr2;

    for (i=n2 -1 ; i>=0 ; i--)
        out_V[i] = 0;

    for (i=n1-1; i >= 0; i--) {
        ptr2 = B[i];
        x = 0;
        for (; ; ptr2++) {
            if ((l = *ptr2) == -1)
                break;
            x ^= in_V[l];
        }
        ptr2 = B[i];
        for (; ; ptr2++) {
            if ((l = *ptr2) == -1)
                break;
            out_V[l] ^= x;
        }
    }

}



static void
VTtimesV(in_1, in_2, out, n2)
unsigned long *in_1;
unsigned long *in_2;
unsigned long *out;
int n2;
{

    register int i,j;
    register unsigned long l,x;
    register unsigned long *in1i, *in2i;
    register int k,m;
    unsigned long mem[WORD_SIZE/BAAT][(1<<BAAT)];

    for (i=0; i < (1<<BAAT); i++)
        for (j=0; j < WORD_SIZE/BAAT; j++)
            mem[j][i] = 0;

    for (i=n2-1,in1i=&(in_1[0]),in2i=&(in_2[0]); i >= 0; i--,in1i++) {
        x = *in2i++;
        if (x) {
            l = *in1i;
            for (j=0; j < WORD_SIZE/BAAT; j++, l>>=BAAT) {
                mem[j][l & ((1<<BAAT)-1)] ^= x;
            }
        }
    }

    for (k=0,m=0; k < WORD_SIZE/BAAT; k++,m+=BAAT) {
        for (l=0,i=1; l <BAAT; l++,i<<=1) {
            x =0;
            for (j=i; j<(1<<BAAT); )
            {
                if ( j & i )
                    x^= mem[k][j++];
                else j+= i;
            }
            out[WORD_SIZEM - m - l] = x;
        }
    }
}



static void
mulNNbyNN(in_1, in_2, out)
unsigned long in_1[WORD_SIZE];
unsigned long in_2[WORD_SIZE];
unsigned long out[WORD_SIZE];
{

    register int i,j;
    register unsigned long bits,m;

    for (i=0; i< WORD_SIZE; i++) {
        bits = 0;
        m = 1UL << WORD_SIZEM;
        for (j=0; j< WORD_SIZE; j++,m >>= 1)
            if (in_1[i] & m)
                bits ^= in_2[j];
        out[i] = bits;
    }

}


static void
addNN(in_1, in_2, out)
unsigned long in_1[WORD_SIZE];
unsigned long in_2[WORD_SIZE];
unsigned long out[WORD_SIZE];

{

    register int i;
    register unsigned long *o, *i1, *i2;

    o = &(out[0]);
    i1 = &(in_1[0]);
    i2 = &(in_2[0]);

    for (i=0; i< WORD_SIZE; i++)
        *o++ = (*i1++) ^ (*i2++);
}



#define swapbits(xx, bi1, bi2) { \
    register unsigned long bit1, bit2,sb1,sb2; \
 \
    sb1 = (1UL<< (WORD_SIZEM - bi1)); \
    sb2 = (1UL<< (WORD_SIZEM - bi2)); \
 \
    bit1 = xx & sb1; \
    bit2 = xx & sb2; \
    xx = (xx & ~(  sb1  | sb2 ) ); \
    if (bit1) xx |= sb2; \
    if (bit2) xx |= sb1; \
}


#define xorbit(xx, bi1, bi2) {     /* bit bi1 ^= bit bi2 */ \
    if (xx &  (1UL<< (WORD_SIZEM - bi2)) ) \
        xx ^= (1UL<<  (WORD_SIZEM - bi1)); \
}



static int
ChooseSi(ViTAVi, Sim1mask, Si, iter)
unsigned long ViTAVi[WORD_SIZE];
unsigned long Sim1mask;
unsigned long Si[WORD_SIZE];
int iter;
{

    unsigned long tmp1[WORD_SIZE];
    int keep[WORD_SIZE+1], ck=0;
    int cols[WORD_SIZE];
    register int i,j,k,ii,x;

    /* mark the columns that must be kept */

    for (i=0; i < WORD_SIZE; i++)
        if ( !(Sim1mask & (1UL << (WORD_SIZEM-i))) ) {
            keep[ck++] = i;
        }

    keep[ck]=0;

    for (i=0; i< WORD_SIZE; i++)
        tmp1[i] = ViTAVi[i];

    for (i=0 ; i < WORD_SIZE; i++)
        cols[i] = i;

    /* do gaussian elimination on cols of tmp1 to find spanning set */

    for (ii=0, i=0; i < WORD_SIZE; i++) {
        /* working on iith column, searching for 1 in the ith row */

        /* first search keep vectors */
        for (j=0; j<ck; j++) {
            /* find column of keep[j] */
            for (k=ii; k<WORD_SIZE; k++)
                if (keep[j]==cols[k])
                    break;
            if (k==WORD_SIZE) {
                printf("buf in ChooseSi - no col %d\n",keep[j]);
                exit(1);
            }
            if (tmp1[i] & (1UL<< (WORD_SIZEM -k)) ) {
                x=k;                    /* got it */
                /* remove keep[j] from list */
                --ck;
                for (k=j; k< ck; k++)
                    keep[k]=keep[k+1];
                j = x;
                goto gotit3;
            }
        }

        if (j==ck) {
            for (j=ii; j < WORD_SIZE; j++)
                if (tmp1[i] & (1UL<< (WORD_SIZEM -j)) )
                    break;
            if (j==WORD_SIZE) {
                continue;       /* no 1 in ith row */
            }
        }

gotit3:
    /* switch columns j and ii */
        if (j != ii) {
            for (k=0; k< WORD_SIZE; k++) {
                swapbits((tmp1[k]), j, ii);
            }
            x= cols[j]; cols[j] = cols[ii]; cols[ii] = x;
        }

        /* now eliminate anything that has a 1 in row i */
        for (j=ii+1; j < WORD_SIZE; j++)
            if (tmp1[i] & (1UL << (WORD_SIZEM -j)) )
                for (k=i; k<WORD_SIZE; k++)
                    xorbit((tmp1[k]), j, ii);

        ii++;
    } /* for (ii=0, i=0; i < WORD_SIZE; i++) */

    if (ck) {
        printf(" ====>  iter %d in ChooseSi()\n",iter);
    }

    /* Si should take cols[0] ... cols[ii-1] */
    for (i=0; i < WORD_SIZE; i++)
        Si[i]=0;

    for (i=0; i < ii; i++) {
        Si[cols[i]] = (1UL << (WORD_SIZEM - i));
    }

    return ii;

}



static void
inverse(in, size, out)                  /* destroys in */
unsigned long in[WORD_SIZE];
int size;
unsigned long out[WORD_SIZE];
{

    register int i,j;
    register unsigned long l, k;

    for (i=0; i < WORD_SIZE; i++)
        out[i] = (1UL << (WORD_SIZEM - i));

    for (i=0; i < size; i++) {
        /* ith row and column */

        /* find a 1 in row j column i */
        l = (1UL << (WORD_SIZEM - i));
        for (j=i; j< size; j++)
            if (in[j] & l )
                break;

        if (j==size) {
            printf("lanczos sub-matrix not invertible: BUG\n");
            exit(1);
        }

        /* swap row j with row i */
        if (j != i) {
            k=in[i]; in[i]= in[j]; in[j]=k;
            k=out[i]; out[i]= out[j]; out[j]=k;
        }

        for (j=i+1; j<size; j++)
            if (in[j] & l) {
                in[j] ^= in[i];
                out[j] ^= out[i];
            }
    }

    for (i=size-2; i >= 0 ; i--)
        for (j=i+1; j<size; j++)
            if (in[i] & (1UL << (WORD_SIZEM - j)))
                out[i] ^= out[j];
}


static void
mulNNbyNNT(in_1, in_2, out)             /* out = in_1 * in_2 transpose */
unsigned long in_1[WORD_SIZE];
unsigned long in_2[WORD_SIZE];
unsigned long out[WORD_SIZE];
{

    unsigned long in_2T[WORD_SIZE];
    register int i,j;
    register unsigned long x,m;

    for (i=0; i < WORD_SIZE; i++)
        in_2T[i] = 0;

    m = (1UL << (WORD_SIZEM ));
    for (i=0; i < WORD_SIZE; i++) {
        x = in_2[i];
        for (j=WORD_SIZEM; j>=0 ; j--) {
            if (x&1)
                in_2T[j] |= m;
            x >>=1;
        }
        m>>=1;
    }

    mulNNbyNN(in_1,in_2T,out);

}


static void
addnN(in_1, in_2, out, n2)
unsigned long *in_1;
unsigned long *in_2;
unsigned long *out;
int n2;
{

    register int i;
    register unsigned long *o, *i1, *i2;

    o = &(out[0]);
    i1 = &(in_1[0]);
    i2 = &(in_2[0]);

    i = n2-1;
    if (n2&1) {
        *o++ = (*i1++) ^ (*i2++);
        i--;
    }

    for (; i>=0; i-=2) {
        *o++ = (*i1++) ^ (*i2++);
        *o++ = (*i1++) ^ (*i2++);
    }

}


static void
mulnNbyNN(in_1, in_2, out, n2)
unsigned long *in_1;
unsigned long *in_2;
unsigned long *out;
int n2;
{

    register int i,k,j;
    register unsigned long bits,l;
    unsigned long mem[(1<<BAAT)*WORD_SIZE/BAAT];


    for (i=0 ;i < WORD_SIZE/BAAT ; i++) {
        l = (i << BAAT); k = i*BAAT;

        mem[l] = 0;
        mem[l|1] = in_2[(BAAT-1)| k];
        mem[l|2] = in_2[(BAAT-2)| k];
        mem[l|3] = in_2[(BAAT-1)| k] ^ in_2[(BAAT-2)|k ];
        mem[l|4] = in_2[(BAAT-3)| k];
        mem[l|5] = in_2[(BAAT-1)| k ] ^ in_2[(BAAT-3)|k ];
        mem[l|6] = in_2[(BAAT-2)| k ] ^ in_2[(BAAT-3)|k ];
        mem[l|7] = in_2[(BAAT-1)| k ] ^ mem[l|6];

        for (j=8; j < (1<<BAAT); ) {
            if (j < 16)
                bits = in_2[(BAAT-4)|k] ^ mem[l|(j - ( 1 << 3))];
            else if (j < 32)
                bits = in_2[(BAAT-5)|k] ^ mem[l|(j - ( 1 << 4))];
            else if (j < 64)
                bits = in_2[(BAAT-6)|k] ^ mem[l|(j - ( 1 << 5))];
            else if (j < 128)
                bits = in_2[(BAAT-7)|k] ^ mem[l|(j - ( 1 << 6))];
            else
                bits = in_2[(BAAT-8)|k] ^ mem[l|(j - ( 1 << 7))];
            mem[l|j++] = bits;
            mem[l|j++] = bits ^ mem[l|1];
            mem[l|j++] = bits ^ mem[l|2];
            mem[l|j++] = bits ^ mem[l|3];
            mem[l|j++] = bits ^ mem[l|4];
            mem[l|j++] = bits ^ mem[l|5];
            mem[l|j++] = bits ^ mem[l|6];
            mem[l|j++] = bits ^ mem[l|7];
        }
    }

    for (i=n2-1; i>=0; i--) {
        l = in_1[i];
        bits = mem[(((WORD_SIZE/BAAT) - 1) <<BAAT) | (l& ((1<<BAAT)-1) )];
        l >>= BAAT;
        switch ( ((WORD_SIZE/BAAT) - 2) ) {
         case 14: 
                bits ^= mem[14*(1<<BAAT) | (l & ((1<<BAAT)-1))]; l>>=BAAT;
                bits ^= mem[13*(1<<BAAT) | (l & ((1<<BAAT)-1))]; l>>=BAAT;
                bits ^= mem[12*(1<<BAAT) | (l & ((1<<BAAT)-1))]; l>>=BAAT;
                bits ^= mem[11*(1<<BAAT) | (l & ((1<<BAAT)-1))]; l>>=BAAT;
                bits ^= mem[10*(1<<BAAT) | (l & ((1<<BAAT)-1))]; l>>=BAAT;
                bits ^= mem[9*(1<<BAAT) | (l & ((1<<BAAT)-1))]; l>>=BAAT;
                bits ^= mem[8*(1<<BAAT) | (l & ((1<<BAAT)-1))]; l>>=BAAT;
                bits ^= mem[7*(1<<BAAT) | (l & ((1<<BAAT)-1))]; l>>=BAAT;
         case 6:
                bits ^= mem[6*(1<<BAAT) | (l & ((1<<BAAT)-1))]; l>>=BAAT;
                bits ^= mem[5*(1<<BAAT) | (l & ((1<<BAAT)-1))]; l>>=BAAT;
                bits ^= mem[4*(1<<BAAT) | (l & ((1<<BAAT)-1))]; l>>=BAAT;
                bits ^= mem[3*(1<<BAAT) | (l & ((1<<BAAT)-1))]; l>>=BAAT;
         case 2:
                bits ^= mem[2*(1<<BAAT) | (l & ((1<<BAAT)-1))]; l>>=BAAT;
                bits ^= mem[1*(1<<BAAT) | (l & ((1<<BAAT)-1))]; l>>=BAAT;
                bits ^= mem[0*(1<<BAAT) | (l & ((1<<BAAT)-1))];
        }

        out[i] = bits;
    }

}


static void
mulNNTbyNN(in_1, in_2, out)
unsigned long in_1[WORD_SIZE];
unsigned long in_2[WORD_SIZE];
unsigned long out[WORD_SIZE];
{

    unsigned long in_1T[WORD_SIZE];
    register int i,j;
    register unsigned long x,m;

    for (i=0; i < WORD_SIZE; i++)
        in_1T[i] = 0;

    m = (1UL << (WORD_SIZEM ));
    for (i=0; i < WORD_SIZE; i++) {
        x = in_1[i];
        for (j=WORD_SIZEM; j>=0 ; j--) {
            if (x&1)
                in_1T[j] |= m;
            x >>=1;
        }
        m>>=1;
    }

    mulNNbyNN(in_1T,in_2,out);

}


#if 0
static int
dimNN(in)
unsigned long in[WORD_SIZE];
{

    register unsigned long k,l;
    register int i,j;

    l = (1UL << (WORD_SIZEM));
    for (i=0; l; i++) {              /* ith row and column */

loop:
        for (j=i; j< WORD_SIZE; j++)
            if (in[j] & l )
                break;

        if (j==WORD_SIZE) {
            l >>=1;
            if (!l) break;
                goto loop;
        }

        /* swap row j with row i */
        if (j != i) {
            k=in[i]; in[i]= in[j]; in[j]=k;
        }

        for (j=i+1; j<WORD_SIZE; j++)
            if (in[j] & l)
                in[j] ^= in[i];

            l >>= 1;
    }

    return  i;

}
#endif

static void
nNcopy(in_V, out_V, n2)
unsigned long in_V[], out_V[];
int n2;
{

    register int j;

    j= n2-1;
    for (; j>=0;j--)
        out_V[j] = in_V[j];

}


#if 0
static void
NNcopy(in, out) 
unsigned long in[], out[];
{

    register int j;

    for (j=WORD_SIZE-1; j>=0;j--)
        out[j] = in[j];

}
#endif

static unsigned long
getmask(diag)
unsigned long diag[];
{

    register int j;
    register unsigned long mask=0;

    for (j=WORD_SIZE-1; j>=0;j--)
        if ((1UL<<j) & diag[WORD_SIZEM-j])
            mask |= (1UL<<j);

    return mask;
}



static void
masknN(in,mask,out, n2)
unsigned long in[], mask, out[];
int n2;
{

    register int j;
    register unsigned long rmask=mask;

    for (j=n2-1; j>=0; j--)
        out[j] = (in[j] & rmask);

}


static void
maskNN(in,mask,out)
unsigned long in[], mask, out[];
{

    register int j;
    register unsigned long rmask=mask;

    for (j=WORD_SIZE-1; j>=0; j--)
        out[j] = (in[j] & rmask);

}


/*
 *  Please see Montgomery's famous block Lanczos paper to understand
 *  variable names and algorithm.
 */
static void
lanczos_main( matrix_fn, size_fb, output_fn )
char *matrix_fn;
char *output_fn;
int size_fb;
{

    int **B;
    int *rel_num;
    int Ni, Nim1 = WORD_SIZE;
    int i, j, k, l;

    unsigned long *X;
    unsigned long *Vi;
    unsigned long *V0;
    unsigned long *AVi;
    unsigned long *Vim1;
    unsigned long *Vim2;
    unsigned long *tmpV;
    unsigned long *tmpV1;

    unsigned long Si[WORD_SIZE];
    unsigned long D[WORD_SIZE];
    unsigned long E[WORD_SIZE];
    unsigned long F[WORD_SIZE];
    unsigned long Winvi[WORD_SIZE];
    unsigned long Winvim1[WORD_SIZE];
    unsigned long Winvim2[WORD_SIZE];
    unsigned long ViTAVi[WORD_SIZE];
    unsigned long ViTAAVi[WORD_SIZE];
    unsigned long SiSiT[WORD_SIZE];
    unsigned long Vim1TAVim1[WORD_SIZE];
    unsigned long tmp1[WORD_SIZE];
    unsigned long tmp2[WORD_SIZE];
    unsigned long tmp3[WORD_SIZE];
    unsigned long Vip1TV0[WORD_SIZE];
    unsigned long ViTV0[WORD_SIZE];
    unsigned long Vim1TV0[WORD_SIZE];
    unsigned long Vim2TV0[WORD_SIZE];
    unsigned long yyy[WORD_SIZE];
    unsigned long xxx[WORD_SIZE];
    unsigned long I_N[WORD_SIZE];

    unsigned long Simask, Sim1mask=0;
    unsigned long xx;
    long seed;
    int vc=0, h;
    FILE *fout, *fseed;
    int n2,n1;
    int sum=0;
    int next;

    /* check sizeof(long) */
    if (SIZEOFLONG != sizeof(long))
      {
	fprintf (stderr, "Error, SIZEOFLONG != sizeof(long)\n");
	exit (1);
      }

    B = ReadB (size_fb, &rel_num, matrix_fn, &n1, &n2 );

    X = (unsigned long *) malloc( n2 * sizeof(unsigned long) );
    Vi = (unsigned long *) malloc( n2 * sizeof(unsigned long) );
    V0 = (unsigned long *) malloc( n2 * sizeof(unsigned long) );
    AVi = (unsigned long *) malloc( n2 * sizeof(unsigned long) );
    Vim1 = (unsigned long *) malloc( n2 * sizeof(unsigned long) );
    Vim2 = (unsigned long *) malloc( n2 * sizeof(unsigned long) );
    tmpV = (unsigned long *) malloc( n2 * sizeof(unsigned long) );
    tmpV1 = (unsigned long *) malloc( n2 * sizeof(unsigned long) );

    fseed = fopen(SEED_FILE,"r");
    if (fseed == NULL)
        seed=12;
    else {
        fscanf(fseed,"%ld",&seed);
        fclose(fseed);
    }
    printf("seed = %ld\n",seed);
    srandom(seed);
    GetY(tmpV,n2);
    for (j=0; j < n2; j++)
        X[j] = tmpV[j];

    /* compute V_0 (i = 0) */
    AtimesV(B,tmpV,Vi,n1,n2);

    for (j=0; j<n2; j++)
        V0[j] = Vi[j];

    VTtimesV(Vi, Vi, ViTV0, n2);

printf("Initializing...\n");    fflush(stdout);

    for (i=0; i<WORD_SIZE; i++)
    /* create  N x N  identity matrix */
        I_N[i] = (1UL << (WORD_SIZEM-i));

    for (i=n2; i ; ) {
        Vim1[--i] = 0, Vim2[i] = 0;
    }

    for (j =0 ; j< WORD_SIZE; j++) {
        Winvim1[j] = 0, Winvim2[j] = 0;
        Sim1mask |= (1UL << j);
        Vim1TAVim1[j] = 0;
        Vim1TV0[j] = 0; Vim2TV0[j] = 0;
    }

    /* used for outputting status: */
    next = n2/WIDTH;

printf("beginning Lanczos iteration\n");
    for (i=0; sum <= n1 ; i++) {
        AtimesV(B,Vi,AVi,n1,n2);

        /* compute (Vi)^T * A * Vi */
        VTtimesV(Vi,AVi, ViTAVi, n2);    /* N by N */

        /* check to see if ViTAVi is 0 */
        for (j=0; j < WORD_SIZE; j++)
            if (ViTAVi[j])
                break;
        if (j==WORD_SIZE)
            break;

        VTtimesV(AVi, AVi, ViTAAVi, n2);

        Ni = ChooseSi(ViTAVi, Sim1mask, Si, i);
        sum += Ni;
if (sum >= next) {
printf(".");
fflush(stdout);
next += (n2/WIDTH);
}

        /* compute Winvi */
        mulNNbyNN(ViTAVi,Si, tmp1);
        mulNNTbyNN(Si,tmp1, tmp2);
        inverse(tmp2, Ni, tmp1);
        mulNNbyNN(Si,tmp1, tmp2);
        mulNNbyNNT(tmp2,Si,Winvi);

        /* Compute D */
        mulNNbyNNT(Si,Si,SiSiT);
        Simask = getmask(SiSiT);
        if (Ni < WORD_SIZE) {
            maskNN(ViTAAVi,Simask,tmp1);
            addNN(tmp1,ViTAVi,yyy);
        }
        else
            addNN(ViTAAVi,ViTAVi,yyy);
        mulNNbyNN(Winvi,yyy,tmp1);
        addNN(I_N, tmp1, D);

        /* Compute E */
        if (Ni < WORD_SIZE) {
            maskNN(ViTAVi,Simask,tmp1);
            mulNNbyNN(Winvim1,tmp1,E);
        }
        else mulNNbyNN(Winvim1,ViTAVi,E);

        /* Compute F */
        if (Nim1 < WORD_SIZE) {
            mulNNbyNN(Vim1TAVim1, Winvim1, tmp1);
            addNN(I_N,tmp1, tmp2);
            mulNNbyNN(Winvim2, tmp2, tmp1);
            if (Ni < WORD_SIZE) {
                mulNNbyNN(tmp1, xxx, tmp3);
                maskNN(tmp3,Simask,F);
            }
            else mulNNbyNN(tmp1, xxx, F);
        }

                /* update X */
        if (i < 3 )
            VTtimesV(Vi,V0,ViTV0, n2);
        mulNNbyNN(Winvi,ViTV0,tmp1);
        mulnNbyNN(Vi,tmp1,tmpV1, n2);

        for (j=n2-1; j >= 0; j--)
            X[j] ^= tmpV1[j];

        /* compute Vip1TV0 */

        if (i > 1) {
            mulNNTbyNN(D,ViTV0, tmp1);
            mulNNTbyNN(E,Vim1TV0,tmp2);
            addNN(tmp1,tmp2,Vip1TV0);
            if (Nim1 < WORD_SIZE) {
                mulNNTbyNN(F,Vim2TV0,tmp1);
                addNN(tmp1,Vip1TV0,Vip1TV0);
            }
        }

        /* Compute the new Vi */

        if (Ni < WORD_SIZE)
            masknN(AVi,Simask,tmpV,n2);
        else
            nNcopy(AVi, tmpV, n2);
        mulnNbyNN(Vi,D,tmpV1, n2);
        addnN(tmpV,tmpV1, tmpV1, n2);
        mulnNbyNN(Vim1,E,tmpV, n2);
        addnN(tmpV,tmpV1,tmpV, n2);
        if (Nim1 < WORD_SIZE) {
            mulnNbyNN(Vim2,F,tmpV1, n2);
            addnN(tmpV,tmpV1,tmpV, n2);
        }
        /* tmpV is the new Vi */

        for (j=n2-1; j>=0; j--) {
            Vim2[j] = Vim1[j];
            Vim1[j] = Vi[j];
            Vi[j] = tmpV[j];
        }

        for (j=0; j< WORD_SIZE; j++) {
            Winvim2[j] = Winvim1[j];
            Winvim1[j] = Winvi[j];
            Vim1TAVim1[j] = ViTAVi[j];
            Vim2TV0[j] = Vim1TV0[j];
            Vim1TV0[j] = ViTV0[j];
            ViTV0[j] = Vip1TV0[j];
            xxx[j] = yyy[j];
        }

        Nim1 = Ni;
        Sim1mask = Simask;
    }


    if (sum > n1) {
        printf("Failed for no apparent reason\n");
        exit(1);
    }

printf("\nEnded at iteration %d (%d orthogonal vectors)\n",i,sum);

    fout = fopen( output_fn, "w" );

    for (j=0; j < n2; j++)
        if (Vi[j])
            goto nonzero2;

#if 0
printf("V is zero\n");
#endif

    /* The vectors in X should be in the null space of B */
    for (i=n1-1; i >= 0; i--) {
        xx = 0;
        for (j=0; ; j++) {
            if ((l = B[i][j]) == -1)
                break;
            xx ^= X[l];
        }
        tmpV[i] = xx;
    }

    /* If the kth column of tmpV is zero, then the kth column of X is in
       the Null Space.  Of course there are much more compact ways to
       output data, but this way it is readable to people.
     */
    for (k=0; k < WORD_SIZE; k++) {
        l = (1UL<< k);
        for (j=0; j<n1; j++)
            if (tmpV[j] & l) break;
        if (j== n1) {
            h=0;
            for (j=0; j<n2; j++)
            if (X[j] & l) {
                fprintf(fout,"%d ",rel_num[j]+1);
                if (++h == 10) {
                    fprintf(fout,"\n");
                    h=0;
                }
            }
            fprintf(fout," -1\n");
            vc++;
        }
    }
    goto free_memory;

nonzero2:

#if 0
printf("V is not zero\n");
#endif

    /* compute B times X and B times Vi */

    for (i=n1-1; i >= 0; i--) {
    xx = 0;
        for (j=0; ; j++) {
            if ((l = B[i][j]) == -1)
                break;
            xx ^= X[l];
        }
        tmpV[i] = xx;
    }

    for (i=n1-1; i >= 0; i--) {
        xx = 0;
        for (j=0; ; j++) {
            if ((l = B[i][j]) == -1)
                break;
            xx ^= Vi[l];
        }
        tmpV1[i] = xx;
    }

    j = 0;
    /* work on ith column of tmpV, X */
    for (i=0; i< WORD_SIZE; i++) {

        /* search for a 1 in row j */
        while (j < n1) {
            for (k=i; k < WORD_SIZE; k++)
                if (tmpV[j] & (1UL << k)) {
                    if (k != i) {
                        for (h=j; h < n1; h++) {
                            nswapbits(tmpV[h],i,k);
                        }
                        for (h=0; h < n2; h++) {
                            nswapbits(X[h],i,k);
                        }
                    }
                    goto gotit;
                }

            for (k=0; k < WORD_SIZE; k++)
                if (tmpV1[j]& (1UL<<k)) {
                    for (h=j; h < n1; h++) {
                        nswapbits2(tmpV1[h],k,tmpV[h],i);
                    }
                    for (h=0; h < n2; h++) {
                        nswapbits2(Vi[h],k, X[h],i);
                    }
                    goto gotit;
                }

            j++;
        }

        if (j==n1) {
            /* found null space vectors */
            for (k=i; k < WORD_SIZE; k++) {
                for (h=0; h <n2; h++)
                    if ((X[h] >> k) & 1)
                        break;
                    if (h < n2) {   /* nonzero */
                        l=0;
                        for (;h <n2; h++)
                            if ((X[h] >> k) & 1) {
                                fprintf(fout,"%d ",rel_num[h]+1);
                                if (++l == 10) {
                                    fprintf(fout,"\n");
                                    l=0;
                                }
                            }
                        fprintf(fout," -1\n");
                        vc++;
                    } /* if (h < n2) */
            } /* for (k=i; k < WORD_SIZE; k++) */
            if (vc) {
                goto done;
            }
        } /* if (j==n1) */

gotit:
    /* eliminate 1's in row j using column i */
        for (k=i+1; k < WORD_SIZE; k++)
            if (tmpV[j] & (1UL<<k)) {
                for (h=j; h < n1; h++)
                    nxorbits(tmpV[h],k,i);
                for (h=0; h < n2; h++)
                    nxorbits(X[h],k,i);
            }
        for (k=0; k < WORD_SIZE; k++)
            if (tmpV1[j] & (1UL<<k)) {
                for (h=j; h < n1; h++)
                    nxorbits2(tmpV1[h],k,tmpV[h],i);
                for (h=0; h < n2; h++)
                    nxorbits2(Vi[h],k,X[h],i);
            }

    } /* for (i=0; i< WORD_SIZE; i++) */

    /* work on ith column of tmpV1, Vi */
    for (i=0; i< WORD_SIZE; i++) {

        /* search for a 1 in row j */
        while (j < n1) {
            for (k=i; k < WORD_SIZE; k++)
                if (tmpV1[j] & (1UL<<k)) {
                    if (k != i) {
                        for (h=j; h < n1; h++) {
                            nswapbits(tmpV1[h],k,i);
                        }
                        for (h=0; h < n2; h++) {
                            nswapbits(Vi[h],k,i);
                        }
                    }
                    goto gotit2;
                }
            j++;
        }

        if (j==n1) {            /* found null space vectors */
            for (k=i; k < WORD_SIZE; k++) {
                for (h=0; h <n2; h++)
                    if (Vi[h] & (1UL<<k))
                        break;
                if (h < n2) {   /* nonzero */
                    l=0;
                    for (;h <n2; h++)
                        if (Vi[h] & (1UL<<k)) {
                            fprintf(fout,"%d ",rel_num[h]+1);
                            if (++l == 10) {
                                fprintf(fout,"\n");
                                l=0;
                            }
                        }
                    fprintf(fout," -1\n");
                    vc++;
                } /* if (h < n2) */
            } /* for (k=i; k < WORD_SIZE; k++) */
            break;
        }

gotit2:
        /* eliminate 1's in row j using column i */
        for (k=i+1; k < WORD_SIZE; k++)
            if (tmpV1[j] & (1UL<<k)) {
                for (h=j; h < n1; h++)
                    nxorbits(tmpV1[h],k,i);
                for (h=0; h < n2; h++)
                    nxorbits(Vi[h],k,i);
            }

    } /* for (i=0; i< WORD_SIZE; i++) */

done:
    printf("output %d vectors\n",vc);


free_memory:
    free( rel_num );
    free( X );
    free( Vi );
    free( V0 );
    free( AVi );
    free( Vim1 );
    free( Vim2 );
    free( tmpV );
    free( tmpV1 );
    for (i=0; i < n1; ++i)
        free( B[i] );
    free( B );

    fclose( fout );
}



static void
usage (void)
{
  fprintf (stderr, "Usage: lanczos -p <params>\n");
  exit (1);
}




/* return non-zero iff a non-trivial factor was found */
static int
process_dependencies (mpz_t N, char *matrix_fn, int *fb, char * deps_fn )
{
  unsigned long i, j;
  int first;
  mpz_t Q; /* collects the large primes */
  mpz_t g;
  int p;
  FILE *fdep, *frel;
  int next_rel_no, curr_rel_no;
  relation rel, r, s;
  int not_done = 1, done;
  int dep_num = 0;

  fdep = fopen( deps_fn , "r");

  /* start reading new dependency */
  while (not_done) {
      frel = fopen (matrix_fn , "r");

      init_relation (r);
      init_relation (s);
      mpz_init_set_ui (Q, 1);
      curr_rel_no = 1;
      first = 1;

      /* get next relation number in dependency: */
      while ( (not_done = fscanf( fdep, "%d", &next_rel_no)) == 1 ) {

          if (next_rel_no == -1)
              /* end of dependency */
              break;

          /* skip past relations not in dependency */
          while (curr_rel_no < next_rel_no) {
              init_relation( rel );
              /* read in relation */
              read_partial (rel, frel, 0);
              clear_relation( rel );
              ++curr_rel_no;
          }

          /* now  cur_rel_no = next_rel_no : process this relation */

          init_relation( rel );
          read_partial (rel, frel, 0);

          if (first) {
              copy_relation (r, rel);
              first = 0;
          }
          else {
              copy_relation (s, r);

              multiply_relation (r, s, rel, N);
          }
          if (rel->q != 0) {
              mpz_mul_ui (Q, Q, rel->q);
              mpz_mod (Q, Q, N);
          }
          

          clear_relation( rel );
          ++curr_rel_no;

      } /* while */

      printf("\nTrying dependency #%d:\n",++dep_num);

      /* check exponents are all even, and take square root */
      for (i = r->n; i-- > 0;) {
          j = r->e[i];
          p = r->p[i];
          if (j % 2) {
              fprintf (stderr, "Error: odd exponent %lu for prime %d\n", j, p);
              exit (1);
              goto end;
          }
          j /= 2;
          while (j--) {
              mpz_mul_si (Q, Q, fb[INDEX(p)]);
              mpz_mod (Q, Q, N);
          }
      }

      mpz_init (g);
      mpz_add (g, r->x, Q);
      mpz_gcd (g, g, N);
      done = (mpz_cmp_ui (g, 1) > 0) && (mpz_cmp (g, N) < 0);
      fprintf (stderr, "gcd(X+Y,N)=");
      mpz_out_str (stderr, 10, g);
      fprintf (stderr, "\n");
      mpz_sub (g, r->x, Q);
      mpz_gcd (g, g, N);
      done = (done || ((mpz_cmp_ui (g, 1)>0) && (mpz_cmp (g, N)<0)));
      fprintf (stderr, "gcd(X-Y,N)=");
      mpz_out_str (stderr, 10, g);
      fprintf (stderr, "\n");
      mpz_clear (g);
      not_done = !done;

end:
      clear_relation (r);
      clear_relation (s);
      mpz_clear (Q);

      fclose( frel );
    }

    fclose( fdep );

  return not_done;
}


static void
factor (N, matrix_fn, fbase_fn, deps_fn )
char *matrix_fn, *fbase_fn, *deps_fn;
mpz_t N;
{
  int *fb;
  unsigned long fbn = 0;

  fbn = count_lines (fbase_fn);
  fb = (int*) malloc ((fbn + 1) * sizeof (int));
  read_factor_base (fb, fbase_fn, fbn);

  process_dependencies( N, matrix_fn, fb, deps_fn );

  free( fb );

}


int
main (int argc, char *argv[])
{
  char *params_fn = NULL, *matrix_fn = NULL, *fbase_fn = NULL;
  char deps_fn[1000];
  int size_fb;
  mpz_t N;
  unsigned long B, LP;

  while (argc > 1)
    {
      if (argc > 2 && strcmp (argv[1], "-p") == 0)
        {
          params_fn = argv[2];
          argv += 2;
          argc -= 2;
        }
      else
        usage ();
    }

  if (params_fn == NULL)
    usage ();

  fbase_fn = (char*) malloc (strlen (params_fn) + strlen (FBASE) + 1);
  strcpy (fbase_fn, params_fn);
  strcat (fbase_fn, FBASE);

  matrix_fn = (char*) malloc (strlen (params_fn) + strlen (MATRIX) + 1);
  strcpy (matrix_fn, params_fn);
  strcat (matrix_fn, MATRIX);

  /* add 1 for "-1" in the factor base */
  size_fb = count_lines (fbase_fn) + 1;

  sprintf( deps_fn, "%s.deps", params_fn );

  lanczos_main( matrix_fn, size_fb, deps_fn );

  printf("%d seconds\n", (int) (cputime()/1000.0 + 0.5));

  mpz_init( N );
  read_params( N, &B, &LP, params_fn );

  factor( N, matrix_fn, fbase_fn, deps_fn );

  printf("total time: %d seconds\n", (int) (cputime()/1000.0 + 0.5));

  free (fbase_fn);
  free (matrix_fn);

  mpz_clear( N );

  return 0;
}
