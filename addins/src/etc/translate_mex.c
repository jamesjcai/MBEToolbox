#include "mex.h"
#include "matrix.h"
#include <stdlib.h>

/* MEX DNA > amino acids
  syntax: [A, Aaa] = translate(str, [icode])
  icode - Genetic code
   1 - Standard;
   2 - Vertebrate Mitochondrial;
   3 - Yeast Mitochondrial;
   4 - Mold, Protozoan, Coelenterate Mitochondrial, and Mycoplasma/Spiroplasma;
   5 - Invertebrate Mitochondrial;
   6 - Ciliate, Dasycladacean, and Hexamita Nuclear */


static const char assign[] = {
	0x0,  0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
	0x0,  0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
	0x0,  0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
	0x0,  0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
	0x0,  0x0, 0x0, 0x1, 0x0, 0x0, 0x0, 0x2, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
	0x0,  0x0, 0x2, 0x0, 0x3, 0x3, 0x0, 0x0, 0x0, 0x3, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
	0x0,  0x0, 0x0, 0x1, 0x0, 0x0, 0x0, 0x2, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
	0x0,  0x0, 0x2, 0x0, 0x3, 0x3, 0x0, 0x0, 0x0, 0x3, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
	0x0,  0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
	0x0,  0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
	0x0,  0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
	0x0,  0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
	0x0,  0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
	0x0,  0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
	0x0,  0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
	0x0,  0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0 };

static const char codon_table[13][65] = {
         {"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVXYXYSSSSXCWCLFLF-"},
         {"KNKNTTTTXSXSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVXYXYSSSSWCWCLFLF-"},
         {"KNKNTTTTRSRSMIMIQHQHPPPPRRRRTTTTEDEDAAAAGGGGVVVVXYXYSSSSWCWCLFLF-"},
         {"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVXYXYSSSSWCWCLFLF-"},
         {"KNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVXYXYSSSSWCWCLFLF-"},
         {"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVQYQYSSSSXCWCLFLF-"},
         {"NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVXYXYSSSSWCWCLFLF-"},
         {"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVXYXYSSSSCCWCLFLF-"},
         {"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVXYXYSSSSXCWCLFLF-"},
         {"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLSLEDEDAAAAGGGGVVVVXYXYSSSSXCWCLFLF-"},
         {"KNKNTTTTGSGSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVXYXYSSSSWCWCLFLF-"},
         {"NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVYYXYSSSSWCWCLFLF-"},
         {"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVXYQYSSSSXCWCLFLF-"}
};


 void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	char	*str    = 0;
	char    *buffer, *bptr;
	char    *buffer2;
	int     c, n, icode = 0;
	int		len1, len2, i;

	if ( nrhs >= 1 && mxIsChar(prhs[0]) ) {

	    len1 = mxGetN(prhs[0]);
		if ((str = mxMalloc(len1 + 1)) == NULL)
			mexErrMsgTxt("translate: mxMalloc(str) error !!!");
		if (mxGetString(prhs[0], str, len1 + 1) != 0)
			mexErrMsgTxt("translate: mxGetString(str) error !!!");	    len2 = len1/3;

		buffer = mxCalloc(len2+1, sizeof(char));
	    bptr   = buffer;

        if (nrhs == 2)
            icode = (int) mxGetScalar(prhs[1]) - 1;
        if ((icode < 0) || (icode > 12))
            mexErrMsgTxt("translate: icode must be an integer from 1 to 13.\n");

		for (i = 0; i < len2; i++) {
		    n = 3 * i;
		    c   = 16*assign[str[n]] + 4*assign[str[n+1]] +assign[str[n+2]];
		    *bptr++ = codon_table[icode][c];
		}

        plhs[0] = mxCreateString(buffer);

        if (nlhs == 2) {
		    buffer2 = mxCalloc(len2*3+1, sizeof(char));
	        bptr   = buffer2;
	        for (i = 0; i < len2; i++) {
	            switch (buffer[i]) {
                case 'A':
                    *bptr++ = 'A'; *bptr++ = 'l'; *bptr++ = 'a'; break;
                case 'C':
                    *bptr++ = 'C'; *bptr++ = 'y'; *bptr++ = 's'; break;
                case 'D':
                    *bptr++ = 'A'; *bptr++ = 's'; *bptr++ = 'p'; break;
                case 'E':
                    *bptr++ = 'G'; *bptr++ = 'l'; *bptr++ = 'u'; break;
                case 'F':
                    *bptr++ = 'P'; *bptr++ = 'h'; *bptr++ = 'e'; break;
                case 'G':
                    *bptr++ = 'G'; *bptr++ = 'l'; *bptr++ = 'y'; break;
                case 'H':
                    *bptr++ = 'H'; *bptr++ = 'i'; *bptr++ = 's'; break;
                case 'I':
                    *bptr++ = 'I'; *bptr++ = 'l'; *bptr++ = 'e'; break;
                case 'K':
                    *bptr++ = 'L'; *bptr++ = 'y'; *bptr++ = 's'; break;
                case 'L':
                    *bptr++ = 'L'; *bptr++ = 'e'; *bptr++ = 'u'; break;
                case 'M':
                    *bptr++ = 'M'; *bptr++ = 'e'; *bptr++ = 't'; break;
                case 'N':
                    *bptr++ = 'A'; *bptr++ = 's'; *bptr++ = 'n'; break;
                case 'P':
                    *bptr++ = 'P'; *bptr++ = 'r'; *bptr++ = 'o'; break;
                case 'Q':
                    *bptr++ = 'G'; *bptr++ = 'l'; *bptr++ = 'n'; break;
                case 'R':
                    *bptr++ = 'A'; *bptr++ = 'r'; *bptr++ = 'g'; break;
                case 'S':
                    *bptr++ = 'S'; *bptr++ = 'e'; *bptr++ = 'r'; break;
                case 'T':
                    *bptr++ = 'T'; *bptr++ = 'h'; *bptr++ = 'r'; break;
                case 'V':
                    *bptr++ = 'V'; *bptr++ = 'a'; *bptr++ = 'l'; break;
                case 'W':
                    *bptr++ = 'T'; *bptr++ = 'r'; *bptr++ = 'p'; break;
                case 'X':
                    *bptr++ = 'A'; *bptr++ = 'm'; *bptr++ = 'b'; break;
                case 'Y':
                    *bptr++ = 'T'; *bptr++ = 'y'; *bptr++ = 'r'; break;
                default:
                    ;
	            }
            }
            plhs[1] = mxCreateString(buffer2);
        }

	} else {
		mexPrintf("usage: [A, Aaa] = translate(cDNA, [icode])");
	}

}
