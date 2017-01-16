/*
 * MATLAB Compiler: 3.0
 * Date: Sun Jun 05 20:26:47 2005
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-x" "-W" "mex" "-L" "C"
 * "-t" "-T" "link:mexlibrary" "libmatlbmx.mlib" "findnest.m" 
 */

#ifndef MLF_V2
#define MLF_V2 1
#endif

#ifndef __findnest_h
#define __findnest_h 1

#ifdef __cplusplus
extern "C" {
#endif

#include "libmatlb.h"

extern void InitializeModule_findnest(void);
extern void TerminateModule_findnest(void);
extern _mexLocalFunctionTable _local_function_table_findnest;

extern mxArray * mlfFindnest(mxArray * * r, mxArray * * s1, mxArray * s);
extern void mlxFindnest(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]);

#ifdef __cplusplus
}
#endif

#endif
