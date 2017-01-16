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

#ifndef __cellstr_h
#define __cellstr_h 1

#ifdef __cplusplus
extern "C" {
#endif

#include "libmatlb.h"

extern void InitializeModule_cellstr(void);
extern void TerminateModule_cellstr(void);
extern _mexLocalFunctionTable _local_function_table_cellstr;

extern mxArray * mlfCellstr(mxArray * s);
extern void mlxCellstr(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]);

#ifdef __cplusplus
}
#endif

#endif
