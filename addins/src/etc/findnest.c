/*
 * MATLAB Compiler: 3.0
 * Date: Sun Jun 05 20:26:47 2005
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-x" "-W" "mex" "-L" "C"
 * "-t" "-T" "link:mexlibrary" "libmatlbmx.mlib" "findnest.m" 
 */
#include "findnest.h"
#include "libmatlbm.h"
static mxArray * _mxarray0_;
static mxArray * _mxarray1_;

static mxChar _array3_[1] = { '(' };
static mxArray * _mxarray2_;

static mxChar _array5_[1] = { ')' };
static mxArray * _mxarray4_;
static mxArray * _mxarray6_;
static mxArray * _mxarray7_;
static mxArray * _mxarray8_;

void InitializeModule_findnest(void) {
    _mxarray0_ = mclInitializeDouble(0.0);
    _mxarray1_ = mclInitializeDouble(1.0);
    _mxarray2_ = mclInitializeString(1, _array3_);
    _mxarray4_ = mclInitializeString(1, _array5_);
    _mxarray6_ = mclInitializeDoubleVector(0, 0, (double *)NULL);
    _mxarray7_ = mclInitializeDouble(2.0);
    _mxarray8_ = mclInitializeDouble(3.0);
}

void TerminateModule_findnest(void) {
    mxDestroyArray(_mxarray8_);
    mxDestroyArray(_mxarray7_);
    mxDestroyArray(_mxarray6_);
    mxDestroyArray(_mxarray4_);
    mxDestroyArray(_mxarray2_);
    mxDestroyArray(_mxarray1_);
    mxDestroyArray(_mxarray0_);
}

static mxArray * Mfindnest(mxArray * * r,
                           mxArray * * s1,
                           int nargout_,
                           mxArray * s);

_mexLocalFunctionTable _local_function_table_findnest
  = { 0, (mexFunctionTableEntry *)NULL };

/*
 * The function "mlfFindnest" contains the normal interface for the "findnest"
 * M-function from file "c:\documents and
 * settings\jamescai\desktop\mbetoolbox\addins\findnest.m" (lines 1-42). This
 * function processes any input arguments and passes them to the implementation
 * version of the function, appearing above.
 */
mxArray * mlfFindnest(mxArray * * r, mxArray * * s1, mxArray * s) {
    int nargout = 1;
    mxArray * l = NULL;
    mxArray * r__ = NULL;
    mxArray * s1__ = NULL;
    mlfEnterNewContext(2, 1, r, s1, s);
    if (r != NULL) {
        ++nargout;
    }
    if (s1 != NULL) {
        ++nargout;
    }
    l = Mfindnest(&r__, &s1__, nargout, s);
    mlfRestorePreviousContext(2, 1, r, s1, s);
    if (r != NULL) {
        mclCopyOutputArg(r, r__);
    } else {
        mxDestroyArray(r__);
    }
    if (s1 != NULL) {
        mclCopyOutputArg(s1, s1__);
    } else {
        mxDestroyArray(s1__);
    }
    return mlfReturnValue(l);
}

/*
 * The function "mlxFindnest" contains the feval interface for the "findnest"
 * M-function from file "c:\documents and
 * settings\jamescai\desktop\mbetoolbox\addins\findnest.m" (lines 1-42). The
 * feval function calls the implementation version of findnest through this
 * function. This function processes any input arguments and passes them to the
 * implementation version of the function, appearing above.
 */
void mlxFindnest(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[1];
    mxArray * mplhs[3];
    int i;
    if (nlhs > 3) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: findnest Line: 1 Column:"
            " 1 The function \"findnest\" was called with m"
            "ore than the declared number of outputs (3)."),
          NULL);
    }
    if (nrhs > 1) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: findnest Line: 1 Column:"
            " 1 The function \"findnest\" was called with m"
            "ore than the declared number of inputs (1)."),
          NULL);
    }
    for (i = 0; i < 3; ++i) {
        mplhs[i] = NULL;
    }
    for (i = 0; i < 1 && i < nrhs; ++i) {
        mprhs[i] = prhs[i];
    }
    for (; i < 1; ++i) {
        mprhs[i] = NULL;
    }
    mlfEnterNewContext(0, 1, mprhs[0]);
    mplhs[0] = Mfindnest(&mplhs[1], &mplhs[2], nlhs, mprhs[0]);
    mlfRestorePreviousContext(0, 1, mprhs[0]);
    plhs[0] = mplhs[0];
    for (i = 1; i < 3 && i < nlhs; ++i) {
        plhs[i] = mplhs[i];
    }
    for (; i < 3; ++i) {
        mxDestroyArray(mplhs[i]);
    }
}

/*
 * The function "Mfindnest" is the implementation version of the "findnest"
 * M-function from file "c:\documents and
 * settings\jamescai\desktop\mbetoolbox\addins\findnest.m" (lines 1-42). It
 * contains the actual compiled code for that M-function. It is a static
 * function and must only be called from one of the interface functions,
 * appearing below.
 */
/*
 * function [l,r,s1] = findnest(s)
 */
static mxArray * Mfindnest(mxArray * * r,
                           mxArray * * s1,
                           int nargout_,
                           mxArray * s) {
    mexLocalFunctionTable save_local_function_table_
      = mclSetCurrentLocalFunctionTable(&_local_function_table_findnest);
    mxArray * l = NULL;
    mxArray * i = NULL;
    mxArray * count = NULL;
    mxArray * ra = NULL;
    mxArray * le = NULL;
    mxArray * address = NULL;
    mxArray * nadr = NULL;
    mxArray * c = NULL;
    mxArray * k = NULL;
    mclCopyArray(&s);
    /*
     * % [l,r,s1] = findnest(s), where s - is a "tree string", 
     * % finds the position of the leftmost terminal cluster
     * % l & r are left and right edges of substring
     * % s1 is the substring found           
     * %
     * % PHYLLAB toolbox v1.0 internal function
     * l = 0; r = 0;
     */
    mlfAssign(&l, _mxarray0_);
    mlfAssign(r, _mxarray0_);
    /*
     * [l,k] = size(s); 
     */
    mlfSize(mlfVarargout(&l, &k, NULL), mclVa(s, "s"), NULL);
    /*
     * c = zeros(1,k); 
     */
    mlfAssign(&c, mlfZeros(_mxarray1_, mclVv(k, "k"), NULL));
    /*
     * nadr = 1; 
     */
    mlfAssign(&nadr, _mxarray1_);
    /*
     * address = zeros(1,k);
     */
    mlfAssign(&address, mlfZeros(_mxarray1_, mclVv(k, "k"), NULL));
    /*
     * % 2, 3 are left and right brackets, respectively
     * le = '('; ra = ')';
     */
    mlfAssign(&le, _mxarray2_);
    mlfAssign(&ra, _mxarray4_);
    /*
     * count = 0;
     */
    mlfAssign(&count, _mxarray0_);
    /*
     * for i = 1:k
     */
    {
        int v_ = mclForIntStart(1);
        int e_ = mclForIntEnd(mclVv(k, "k"));
        if (v_ > e_) {
            mlfAssign(&i, _mxarray6_);
        } else {
            /*
             * if     strcmp(s(i),le) == 1
             * c(nadr) = 2;
             * address(nadr)=i;
             * nadr = nadr + 1;
             * count = count+1;
             * elseif strcmp(s(i),ra) == 1
             * c(nadr) = 3; 
             * address(nadr)=i;
             * nadr = nadr + 1;
             * count = count+1;
             * end 
             * end 
             */
            for (; ; ) {
                if (mclEqBool(
                      mlfStrcmp(
                        mclIntArrayRef1(mclVa(s, "s"), v_), mclVv(le, "le")),
                      _mxarray1_)) {
                    mclArrayAssign1(&c, _mxarray7_, mclVv(nadr, "nadr"));
                    mclArrayAssign1(
                      &address, mlfScalar(v_), mclVv(nadr, "nadr"));
                    mlfAssign(&nadr, mclPlus(mclVv(nadr, "nadr"), _mxarray1_));
                    mlfAssign(
                      &count, mclPlus(mclVv(count, "count"), _mxarray1_));
                } else if (mclEqBool(
                             mlfStrcmp(
                               mclIntArrayRef1(mclVa(s, "s"), v_),
                               mclVv(ra, "ra")),
                             _mxarray1_)) {
                    mclArrayAssign1(&c, _mxarray8_, mclVv(nadr, "nadr"));
                    mclArrayAssign1(
                      &address, mlfScalar(v_), mclVv(nadr, "nadr"));
                    mlfAssign(&nadr, mclPlus(mclVv(nadr, "nadr"), _mxarray1_));
                    mlfAssign(
                      &count, mclPlus(mclVv(count, "count"), _mxarray1_));
                }
                if (v_ == e_) {
                    break;
                }
                ++v_;
            }
            mlfAssign(&i, mlfScalar(v_));
        }
    }
    /*
     * if count < 2
     */
    if (mclLtBool(mclVv(count, "count"), _mxarray7_)) {
        /*
         * r=0; %%%%%%%%%%%%%%%%%%
         */
        mlfAssign(r, _mxarray0_);
        /*
         * s1 = s; %%%%%%%%%%%%%%%
         */
        mlfAssign(s1, mclVa(s, "s"));
        /*
         * return;
         */
        goto return_;
    /*
     * end
     */
    }
    /*
     * for i = 2:nadr
     */
    {
        int v_ = mclForIntStart(2);
        int e_ = mclForIntEnd(mclVv(nadr, "nadr"));
        if (v_ > e_) {
            mlfAssign(&i, _mxarray6_);
        } else {
            /*
             * if  c(i-1) == 2 & c(i) == 3
             * l = address(i-1);
             * r = address(i);
             * break;
             * end
             * end
             */
            for (; ; ) {
                mxArray * a_
                  = mclInitialize(
                      mclEq(
                        mclIntArrayRef1(mclVv(c, "c"), v_ - 1), _mxarray7_));
                if (mlfTobool(a_)
                    && mlfTobool(
                         mclAnd(
                           a_,
                           mclEq(
                             mclIntArrayRef1(mclVv(c, "c"), v_),
                             _mxarray8_)))) {
                    mxDestroyArray(a_);
                    mlfAssign(
                      &l, mclIntArrayRef1(mclVv(address, "address"), v_ - 1));
                    mlfAssign(
                      r, mclIntArrayRef1(mclVv(address, "address"), v_));
                    break;
                } else {
                    mxDestroyArray(a_);
                }
                if (v_ == e_) {
                    break;
                }
                ++v_;
            }
            mlfAssign(&i, mlfScalar(v_));
        }
    }
    /*
     * address = []; c = [];        
     */
    mlfAssign(&address, _mxarray6_);
    mlfAssign(&c, _mxarray6_);
    /*
     * s1 = s(l:r);
     */
    mlfAssign(
      s1,
      mclArrayRef1(
        mclVa(s, "s"), mlfColon(mclVv(l, "l"), mclVv(*r, "r"), NULL)));
    return_:
    mclValidateOutput(l, 1, nargout_, "l", "findnest");
    mclValidateOutput(*r, 2, nargout_, "r", "findnest");
    mclValidateOutput(*s1, 3, nargout_, "s1", "findnest");
    mxDestroyArray(k);
    mxDestroyArray(c);
    mxDestroyArray(nadr);
    mxDestroyArray(address);
    mxDestroyArray(le);
    mxDestroyArray(ra);
    mxDestroyArray(count);
    mxDestroyArray(i);
    mxDestroyArray(s);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return l;
}
