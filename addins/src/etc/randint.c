/* Function: RandInt() *************************************************
 *

 * Purpose:  Generates random integer matrix
 *
 * syntax R = randint(min,
max, m, n, [prior])          
 * Args:     min   - minimum integer
 *
max   - maximum integer
 *           m, n  - generates (m x n) matrix
 *
prior - prior frequencies of the integers [1x(max-min+1)]
 * version
1.0, 19:26, 25-Aug-2004, created by m.tada

***********************************************************************/

#include <stdlib.h>
#include <time.h>
#include "mex.h"
 
/* MAIN FUNCTION
*/
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray
*prhs[] ) {

    int     m, n;                  /* output matrix size
*/
    int     minN, maxN, N = 0;     /* minimum and maximum number
*/
    int     *table, *pt, size = 0; /* table size and its size      */

int     i, j, val;             /* generic counter              */

double  *mat;                  /* pointer to the output matrix */
    double
*freq, *fptr, sum = 0;
    int     *ifreq; 
    

	/* Check for
proper number of input arguments. */
	if ( nrhs < 2 ) {

	mexErrMsgTxt("Too few input arguments");
	} else if (
nrhs > 5 ) {
		mexErrMsgTxt("Too many input arguments");
	}

	
	if ( nlhs > 1 )
		mexErrMsgTxt("Too many output arguments");
    
	minN = (int) mxGetScalar(prhs[0]);	

	maxN = (int) mxGetScalar(prhs[1]);
	
	if ( nrhs == 2 ) {

	    m = 1;
	    n = 1;
	} else if (nrhs == 3) {

m = n = mxGetScalar(prhs[2]);
	} else {
	    m = mxGetScalar(prhs[2]);
	    n = mxGetScalar(prhs[3]);
	}
    
    if (nrhs== 5) {
        if ((N = mxGetNumberOfElements(prhs[4])) == maxN-minN+1)  {

fptr  = mxGetPr(prhs[4]);
            freq  = mxMalloc( N * sizeof(double) );
            ifreq = mxMalloc( N * sizeof(int) );

for (i = 0; i < N; i++)
                sum += freq[i] = *fptr++;

for (i = 0; i < N; i++)
                freq[i] /= sum;
            for (i =0; i < N; i++) {
                size += ifreq[i] = (int) 1000*freq[i];

}
            table = mxMalloc( size * sizeof(int) );

pt    = table;
            for (i = 0; i < N; i++) {
                for (j
= 0; j < ifreq[i]; j++)
                    *pt++ = minN + i;
            }

} else {
            mexPrintf("Length of frequency vector should be %d...\nUniform distribution is used.\n", maxN-minN+1);
            N = 0;

}
    }
    
	plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
	mat     =
mxGetPr(plhs[0]);
	
	srand( time(NULL) + rand() );
	
    if
(N) {
        for (i = 0; i < (m * n); i++) {
            j      = (int) (
(rand() * size) / (RAND_MAX + 1.0) );
            val    = table[j];

*mat++ = (double) val;
        }
    } else {
        for (i = 0; i < (m *
n); i++) {
            val    = (int) ( (maxN - minN + 1.0)*( 1.0 /
(RAND_MAX + 1.0) ) 
                              * rand()  + minN );

*mat++ = (double) val;
        }
    }   
}