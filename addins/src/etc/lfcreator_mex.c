/*
 * MATLAB Compiler: 3.0
 * Date: Sun Jun 05 20:29:16 2005
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-x" "-W" "mex" "-L" "C"
 * "-t" "-T" "link:mexlibrary" "libmatlbmx.mlib" "lfcreator.m" 
 */

#ifndef MLF_V2
#define MLF_V2 1
#endif

#include "libmatlb.h"
#include "lfcreator.h"
#include "cellstr.h"
#include "num2str.h"

extern _mex_information _mex_info;

static mexFunctionTableEntry function_table[1]
  = { { "lfcreator", mlxLfcreator, 2, 1, &_local_function_table_lfcreator } };

static _mexInitTermTableEntry init_term_table[1]
  = { { InitializeModule_lfcreator, TerminateModule_lfcreator } };

/*
 * The function "Mnum2str" is the MATLAB callback version of the "num2str"
 * function from file "d:\matlab6p5\toolbox\matlab\strfun\num2str.m". It
 * performs a callback to MATLAB to run the "num2str" function, and passes any
 * resulting output arguments back to its calling function.
 */
static mxArray * Mnum2str(int nargout_, mxArray * x, mxArray * f) {
    mxArray * s = NULL;
    mclFevalCallMATLAB(
      mclNVarargout(nargout_, 0, &s, NULL), "num2str", x, f, NULL);
    return s;
}

/*
 * The function "Mcellstr" is the MATLAB callback version of the "cellstr"
 * function from file "d:\matlab6p5\toolbox\matlab\strfun\cellstr.m". It
 * performs a callback to MATLAB to run the "cellstr" function, and passes any
 * resulting output arguments back to its calling function.
 */
static mxArray * Mcellstr(int nargout_, mxArray * s) {
    mxArray * c = NULL;
    mclFevalCallMATLAB(
      mclNVarargout(nargout_, 0, &c, NULL), "cellstr", s, NULL);
    return c;
}

/*
 * The function "mexLibrary" is a Compiler-generated mex wrapper, suitable for
 * building a MEX-function. It initializes any persistent variables as well as
 * a function table for use by the feval function. It then calls the function
 * "mlxLfcreator". Finally, it clears the feval table and exits.
 */
mex_information mexLibrary(void) {
    mclMexLibraryInit();
    return &_mex_info;
}

_mex_information _mex_info
  = { 1, 1, function_table, 0, NULL, 0, NULL, 1, init_term_table };

/*
 * The function "mlfNum2str" contains the normal interface for the "num2str"
 * M-function from file "d:\matlab6p5\toolbox\matlab\strfun\num2str.m" (lines
 * 0-0). This function processes any input arguments and passes them to the
 * implementation version of the function, appearing above.
 */
mxArray * mlfNum2str(mxArray * x, mxArray * f) {
    int nargout = 1;
    mxArray * s = NULL;
    mlfEnterNewContext(0, 2, x, f);
    s = Mnum2str(nargout, x, f);
    mlfRestorePreviousContext(0, 2, x, f);
    return mlfReturnValue(s);
}

/*
 * The function "mlxNum2str" contains the feval interface for the "num2str"
 * M-function from file "d:\matlab6p5\toolbox\matlab\strfun\num2str.m" (lines
 * 0-0). The feval function calls the implementation version of num2str through
 * this function. This function processes any input arguments and passes them
 * to the implementation version of the function, appearing above.
 */
void mlxNum2str(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[2];
    mxArray * mplhs[1];
    int i;
    if (nlhs > 1) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: num2str Line: 1 Column: "
            "1 The function \"num2str\" was called with mor"
            "e than the declared number of outputs (1)."),
          NULL);
    }
    if (nrhs > 2) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: num2str Line: 1 Column:"
            " 1 The function \"num2str\" was called with m"
            "ore than the declared number of inputs (2)."),
          NULL);
    }
    for (i = 0; i < 1; ++i) {
        mplhs[i] = NULL;
    }
    for (i = 0; i < 2 && i < nrhs; ++i) {
        mprhs[i] = prhs[i];
    }
    for (; i < 2; ++i) {
        mprhs[i] = NULL;
    }
    mlfEnterNewContext(0, 2, mprhs[0], mprhs[1]);
    mplhs[0] = Mnum2str(nlhs, mprhs[0], mprhs[1]);
    mlfRestorePreviousContext(0, 2, mprhs[0], mprhs[1]);
    plhs[0] = mplhs[0];
}

/*
 * The function "mlfCellstr" contains the normal interface for the "cellstr"
 * M-function from file "d:\matlab6p5\toolbox\matlab\strfun\cellstr.m" (lines
 * 0-0). This function processes any input arguments and passes them to the
 * implementation version of the function, appearing above.
 */
mxArray * mlfCellstr(mxArray * s) {
    int nargout = 1;
    mxArray * c = NULL;
    mlfEnterNewContext(0, 1, s);
    c = Mcellstr(nargout, s);
    mlfRestorePreviousContext(0, 1, s);
    return mlfReturnValue(c);
}

/*
 * The function "mlxCellstr" contains the feval interface for the "cellstr"
 * M-function from file "d:\matlab6p5\toolbox\matlab\strfun\cellstr.m" (lines
 * 0-0). The feval function calls the implementation version of cellstr through
 * this function. This function processes any input arguments and passes them
 * to the implementation version of the function, appearing above.
 */
void mlxCellstr(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[1];
    mxArray * mplhs[1];
    int i;
    if (nlhs > 1) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: cellstr Line: 1 Column: "
            "1 The function \"cellstr\" was called with mor"
            "e than the declared number of outputs (1)."),
          NULL);
    }
    if (nrhs > 1) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: cellstr Line: 1 Column:"
            " 1 The function \"cellstr\" was called with m"
            "ore than the declared number of inputs (1)."),
          NULL);
    }
    for (i = 0; i < 1; ++i) {
        mplhs[i] = NULL;
    }
    for (i = 0; i < 1 && i < nrhs; ++i) {
        mprhs[i] = prhs[i];
    }
    for (; i < 1; ++i) {
        mprhs[i] = NULL;
    }
    mlfEnterNewContext(0, 1, mprhs[0]);
    mplhs[0] = Mcellstr(nlhs, mprhs[0]);
    mlfRestorePreviousContext(0, 1, mprhs[0]);
    plhs[0] = mplhs[0];
}
