/* --------------------------------------------------- */
/* Automatically generated by mwrap                    */
/* --------------------------------------------------- */

/* Code generated by mwrap */
/*
  Copyright statement for mwrap:

  mwrap -- MEX file generation for MATLAB and Octave
  Copyright (c) 2007-2008 David Bindel

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in
  all copies or substantial portions of the Software.

  You may distribute a work that contains part or all of the source code
  generated by mwrap under the terms of your choice.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
  THE SOFTWARE.
*/

#include <mex.h>
#include <stdio.h>
#include <string.h>


#ifndef ulong
#  define ulong unsigned long
#endif
#ifndef uint
#  define uint  unsigned int
#endif
#ifndef uchar
#  define uchar unsigned char
#endif


/*
 * Support for 32-bit and 64-bit MEX files
 */
#ifndef mwSize
#  define mwSize int
#endif
#ifndef mwIndex
#  define mwIndex int
#endif
#ifndef mwSignedIndex
#  define mwSignedIndex int
#endif


/*
 * Records for call profile.
 */
int* mexprofrecord_= NULL;


/*
 * Support routines for copying data into and out of the MEX stubs
 */

void* mxWrapGetP(const mxArray* a, const char* fmt, const char** e)
{
    void* p = 0;
    mxArray* ap;
    if (mxGetClassID(a) == mxDOUBLE_CLASS && 
        mxGetM(a)*mxGetN(a) == 1 && *mxGetPr(a) == 0)
        return p;
    if (mxIsChar(a)) {
        char pbuf[128];
        mxGetString(a, pbuf, sizeof(pbuf));
        sscanf(pbuf, fmt, &p);
    } 
#ifdef R2008OO
    else if (ap = mxGetProperty(a, 0, "mwptr")) {
        return mxWrapGetP(ap, fmt, e);
    }
#endif
    if (p == 0)
        *e = "Invalid pointer";
    return p;
}

mxArray* mxWrapCreateP(void* p, const char* fmt)
{
    if (p == 0) {
        mxArray* z = mxCreateDoubleMatrix(1,1, mxREAL);
        *mxGetPr(z) = 0;
        return z;
    } else {
        char pbuf[128];
        sprintf(pbuf, fmt, p);
        return mxCreateString(pbuf);
    }
}

mxArray* mxWrapStrncpy(const char* s)
{
    if (s) {
        return mxCreateString(s);
    } else {
        mxArray* z = mxCreateDoubleMatrix(1,1, mxREAL);
        *mxGetPr(z) = 0;
        return z;
    }
}

double mxWrapGetScalar(const mxArray* a, const char** e)
{
    if (!a || mxGetClassID(a) != mxDOUBLE_CLASS || mxGetM(a)*mxGetN(a) != 1) {
        *e = "Invalid scalar argument";
        return 0;
    }
    return *mxGetPr(a);
}

char* mxWrapGetString(const mxArray* a, const char** e)
{
    char* s;
    mwSize slen;
    if (!a || (!mxIsChar(a) && mxGetM(a)*mxGetN(a) > 0)) {
        *e = "Invalid string argument";
        return NULL;
    }
    slen = mxGetM(a)*mxGetN(a) + 1;
    s = (char*) mxMalloc(slen);
    if (mxGetM(a)*mxGetN(a) == 0)
        *s = 0;
    else
        mxGetString(a, s, slen);
    return s;
}


#define mxWrapGetArrayDef(func, T) \
T* func(const mxArray* a, const char** e)     \
{ \
    T* array; \
    mwSize arraylen; \
    mwIndex i; \
    T* p; \
    double* q; \
    if (!a || mxGetClassID(a) != mxDOUBLE_CLASS) { \
        *e = "Invalid array argument"; \
        return 0; \
    } \
    arraylen = mxGetM(a)*mxGetN(a); \
    array = (T*) mxMalloc(mxGetM(a)*mxGetN(a) * sizeof(T)); \
    p = array; \
    q = mxGetPr(a); \
    for (i = 0; i < arraylen; ++i) \
        *p++ = (T) (*q++); \
    return array; \
}


#define mxWrapCopyDef(func, T) \
void func(mxArray* a, const T* q, mwSize n) \
{ \
    mwIndex i; \
    double* p = mxGetPr(a); \
    for (i = 0; i < n; ++i) \
        *p++ = *q++; \
}


#define mxWrapReturnDef(func, T) \
mxArray* func(const T* q, mwSize m, mwSize n) \
{ \
    mwIndex i; \
    double* p; \
    if (!q) { \
        return mxCreateDoubleMatrix(0,0, mxREAL); \
    } else { \
        mxArray* a = mxCreateDoubleMatrix(m,n, mxREAL); \
        p = mxGetPr(a); \
        for (i = 0; i < m*n; ++i) \
            *p++ = *q++; \
        return a; \
    } \
}


#define mxWrapGetScalarZDef(func, T, ZT, setz) \
void func(T* z, const mxArray* a) \
{ \
    double* pr = mxGetPr(a); \
    double* pi = mxGetPi(a); \
    setz(z, (ZT) *pr, (pi ? (ZT) *pi : (ZT) 0)); \
}


#define mxWrapGetArrayZDef(func, T, ZT, setz) \
T* func(const mxArray* a, const char** e) \
{ \
    T* array; \
    mwSize arraylen; \
    mwIndex i; \
    T* p; \
    double* qr; \
    double* qi; \
    if (!a || mxGetClassID(a) != mxDOUBLE_CLASS) { \
        *e = "Invalid array argument"; \
        return 0; \
    } \
    arraylen = mxGetM(a)*mxGetN(a); \
    array = (T*) mxMalloc(mxGetM(a)*mxGetN(a) * sizeof(T)); \
    p = array; \
    qr = mxGetPr(a); \
    qi = mxGetPi(a); \
    for (i = 0; i < arraylen; ++i) { \
        ZT val_qr = *qr++; \
        ZT val_qi = (qi ? (ZT) *qi++ : (ZT) 0); \
        setz(p, val_qr, val_qi); \
        ++p; \
    } \
    return array; \
}


#define mxWrapCopyZDef(func, T, real, imag) \
void func(mxArray* a, const T* q, mwSize n) \
{ \
    mwIndex i; \
    double* pr = mxGetPr(a); \
    double* pi = mxGetPi(a); \
    for (i = 0; i < n; ++i) { \
        *pr++ = real(*q); \
        *pi++ = imag(*q); \
        ++q; \
    } \
}


#define mxWrapReturnZDef(func, T, real, imag) \
mxArray* func(const T* q, mwSize m, mwSize n) \
{ \
    mwIndex i; \
    double* pr; \
    double* pi; \
    if (!q) { \
        return mxCreateDoubleMatrix(0,0, mxCOMPLEX); \
    } else { \
        mxArray* a = mxCreateDoubleMatrix(m,n, mxCOMPLEX); \
        pr = mxGetPr(a); \
        pi = mxGetPi(a); \
        for (i = 0; i < m*n; ++i) { \
            *pr++ = real(*q); \
            *pi++ = imag(*q); \
            ++q; \
        } \
        return a; \
    } \
}

#include <complex.h>

typedef _Complex double dcomplex;
#define real_dcomplex(z) creal(z)
#define imag_dcomplex(z) cimag(z)
#define setz_dcomplex(z,r,i)  *z = r + i*_Complex_I

typedef _Complex float fcomplex;
#define real_fcomplex(z) crealf(z)
#define imag_fcomplex(z) cimagf(z)
#define setz_fcomplex(z,r,i)  *z = r + i*_Complex_I

/* Array copier definitions */
mxWrapGetArrayDef(mxWrapGetArray_bool, bool)
mxWrapCopyDef    (mxWrapCopy_bool,     bool)
mxWrapReturnDef  (mxWrapReturn_bool,   bool)
mxWrapGetArrayDef(mxWrapGetArray_char, char)
mxWrapCopyDef    (mxWrapCopy_char,     char)
mxWrapReturnDef  (mxWrapReturn_char,   char)
mxWrapGetArrayDef(mxWrapGetArray_double, double)
mxWrapCopyDef    (mxWrapCopy_double,     double)
mxWrapReturnDef  (mxWrapReturn_double,   double)
mxWrapGetArrayDef(mxWrapGetArray_float, float)
mxWrapCopyDef    (mxWrapCopy_float,     float)
mxWrapReturnDef  (mxWrapReturn_float,   float)
mxWrapGetArrayDef(mxWrapGetArray_int, int)
mxWrapCopyDef    (mxWrapCopy_int,     int)
mxWrapReturnDef  (mxWrapReturn_int,   int)
mxWrapGetArrayDef(mxWrapGetArray_long, long)
mxWrapCopyDef    (mxWrapCopy_long,     long)
mxWrapReturnDef  (mxWrapReturn_long,   long)
mxWrapGetArrayDef(mxWrapGetArray_mwIndex, mwIndex)
mxWrapCopyDef    (mxWrapCopy_mwIndex,     mwIndex)
mxWrapReturnDef  (mxWrapReturn_mwIndex,   mwIndex)
mxWrapGetArrayDef(mxWrapGetArray_mwSignedIndex, mwSignedIndex)
mxWrapCopyDef    (mxWrapCopy_mwSignedIndex,     mwSignedIndex)
mxWrapReturnDef  (mxWrapReturn_mwSignedIndex,   mwSignedIndex)
mxWrapGetArrayDef(mxWrapGetArray_mwSize, mwSize)
mxWrapCopyDef    (mxWrapCopy_mwSize,     mwSize)
mxWrapReturnDef  (mxWrapReturn_mwSize,   mwSize)
mxWrapGetArrayDef(mxWrapGetArray_size_t, size_t)
mxWrapCopyDef    (mxWrapCopy_size_t,     size_t)
mxWrapReturnDef  (mxWrapReturn_size_t,   size_t)
mxWrapGetArrayDef(mxWrapGetArray_uchar, uchar)
mxWrapCopyDef    (mxWrapCopy_uchar,     uchar)
mxWrapReturnDef  (mxWrapReturn_uchar,   uchar)
mxWrapGetArrayDef(mxWrapGetArray_uint, uint)
mxWrapCopyDef    (mxWrapCopy_uint,     uint)
mxWrapReturnDef  (mxWrapReturn_uint,   uint)
mxWrapGetArrayDef(mxWrapGetArray_ulong, ulong)
mxWrapCopyDef    (mxWrapCopy_ulong,     ulong)
mxWrapReturnDef  (mxWrapReturn_ulong,   ulong)
mxWrapGetScalarZDef(mxWrapGetScalar_fcomplex, fcomplex,
                    float, setz_fcomplex)
mxWrapGetArrayZDef (mxWrapGetArray_fcomplex, fcomplex,
                    float, setz_fcomplex)
mxWrapCopyZDef     (mxWrapCopy_fcomplex, fcomplex,
                    real_fcomplex, imag_fcomplex)
mxWrapReturnZDef   (mxWrapReturn_fcomplex, fcomplex,
                    real_fcomplex, imag_fcomplex)
mxWrapGetScalarZDef(mxWrapGetScalar_dcomplex, dcomplex,
                    double, setz_dcomplex)
mxWrapGetArrayZDef (mxWrapGetArray_dcomplex, dcomplex,
                    double, setz_dcomplex)
mxWrapCopyZDef     (mxWrapCopy_dcomplex, dcomplex,
                    real_dcomplex, imag_dcomplex)
mxWrapReturnZDef   (mxWrapReturn_dcomplex, dcomplex,
                    real_dcomplex, imag_dcomplex)

#if defined(MWF77_CAPS)
#define MWF77_interpmatnoalloc INTERPMATNOALLOC
#define MWF77_extrapnoalloc EXTRAPNOALLOC
#elif defined(MWF77_UNDERSCORE1)
#define MWF77_interpmatnoalloc __dspline_MOD_interpmatnoalloc
#define MWF77_extrapnoalloc __dspline_MOD_extrapnoalloc
#else /* f2c convention */
#define MWF77_interpmatnoalloc __dspline_MOD_interpmatnoalloc
#define MWF77_extrapnoalloc __dspline_MOD_extrapnoalloc
#endif

#ifdef __cplusplus
extern "C" { /* Prevent C++ name mangling */
#endif

#ifndef MWF77_RETURN
#define MWF77_RETURN int
#endif

MWF77_RETURN MWF77_interpmatnoalloc(int*, double*, double*, int*, int*, int*, double*, double*, int*);
MWF77_RETURN MWF77_extrapnoalloc(double*, int*);

#ifdef __cplusplus
} /* end extern C */
#endif

/* ---- gateway.mw: 30 ----
 * interpmatnoalloc(int r, double[] tinterp, double dt, int m, output int[1] jmax, output int[1] jmin, inout double[] umat, inout double[] upmat, int maxinds);
 */
const char* stubids1_ = "interpmatnoalloc(i int, i double[], i double, i int, o int[x], o int[x], io double[], io double[], i int)";

void mexStub1(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    int         in0_;    /* r          */
    double*     in1_ =0; /* tinterp    */
    double      in2_;    /* dt         */
    int         in3_;    /* m          */
    double*     in4_ =0; /* umat       */
    double*     in5_ =0; /* upmat      */
    int         in6_;    /* maxinds    */
    int*        out0_=0; /* jmax       */
    int*        out1_=0; /* jmin       */
    mwSize      dim7_;   /* 1          */
    mwSize      dim8_;   /* 1          */

    dim7_ = (mwSize) mxWrapGetScalar(prhs[7], &mw_err_txt_);
    dim8_ = (mwSize) mxWrapGetScalar(prhs[8], &mw_err_txt_);

    in0_ = (int) mxWrapGetScalar(prhs[0], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[1])*mxGetN(prhs[1]) != 0) {
        in1_ = mxGetPr(prhs[1]);
    } else
        in1_ = NULL;
    in2_ = (double) mxWrapGetScalar(prhs[2], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    in3_ = (int) mxWrapGetScalar(prhs[3], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != 0) {
        in4_ = mxWrapGetArray_double(prhs[4], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in4_ = NULL;
    if (mxGetM(prhs[5])*mxGetN(prhs[5]) != 0) {
        in5_ = mxWrapGetArray_double(prhs[5], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in5_ = NULL;
    in6_ = (int) mxWrapGetScalar(prhs[6], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    out0_ = (int*) mxMalloc(dim7_*sizeof(int));
    out1_ = (int*) mxMalloc(dim8_*sizeof(int));
    if (mexprofrecord_)
        mexprofrecord_[1]++;
    MWF77_interpmatnoalloc(&in0_, in1_, &in2_, &in3_, out0_, out1_, in4_, in5_, &in6_);
    plhs[0] = mxCreateDoubleMatrix(dim7_, 1, mxREAL);
    mxWrapCopy_int(plhs[0], out0_, dim7_);
    plhs[1] = mxCreateDoubleMatrix(dim8_, 1, mxREAL);
    mxWrapCopy_int(plhs[1], out1_, dim8_);
    plhs[2] = mxCreateDoubleMatrix(mxGetM(prhs[4]), mxGetN(prhs[4]), mxREAL);
    mxWrapCopy_double(plhs[2], in4_, mxGetM(prhs[4])*mxGetN(prhs[4]));
    plhs[3] = mxCreateDoubleMatrix(mxGetM(prhs[5]), mxGetN(prhs[5]), mxREAL);
    mxWrapCopy_double(plhs[3], in5_, mxGetM(prhs[5])*mxGetN(prhs[5]));

mw_err_label:
    if (out0_) mxFree(out0_);
    if (out1_) mxFree(out1_);
    if (in4_)  mxFree(in4_);
    if (in5_)  mxFree(in5_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ---- gateway.mw: 53 ----
 * extrapnoalloc(inout double[] xcof, int m);
 */
const char* stubids2_ = "extrapnoalloc(io double[], i int)";

void mexStub2(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    double*     in0_ =0; /* xcof       */
    int         in1_;    /* m          */

    if (mxGetM(prhs[0])*mxGetN(prhs[0]) != 0) {
        in0_ = mxWrapGetArray_double(prhs[0], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in0_ = NULL;
    in1_ = (int) mxWrapGetScalar(prhs[1], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mexprofrecord_)
        mexprofrecord_[2]++;
    MWF77_extrapnoalloc(in0_, &in1_);
    plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[0]), mxGetN(prhs[0]), mxREAL);
    mxWrapCopy_double(plhs[0], in0_, mxGetM(prhs[0])*mxGetN(prhs[0]));

mw_err_label:
    if (in0_)  mxFree(in0_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ----
 */
void mexFunction(int nlhs, mxArray* plhs[],
                 int nrhs, const mxArray* prhs[])
{
    char id[512];
    if (nrhs == 0) {
        mexPrintf("Mex function installed\n");
        return;
    }

    if (mxGetString(prhs[0], id, sizeof(id)) != 0)
        mexErrMsgTxt("Identifier should be a string");
    else if (strcmp(id, stubids1_) == 0)
        mexStub1(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, stubids2_) == 0)
        mexStub2(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, "*profile on*") == 0) {
        if (!mexprofrecord_) {
            mexprofrecord_ = (int*) malloc(3 * sizeof(int));
            mexLock();
        }
        memset(mexprofrecord_, 0, 3 * sizeof(int));
    } else if (strcmp(id, "*profile off*") == 0) {
        if (mexprofrecord_) {
            free(mexprofrecord_);
            mexUnlock();
        }
        mexprofrecord_ = NULL;
    } else if (strcmp(id, "*profile report*") == 0) {
        if (!mexprofrecord_)
            mexPrintf("Profiler inactive\n");
        mexPrintf("%d calls to gateway.mw:30\n", mexprofrecord_[1]);
        mexPrintf("%d calls to gateway.mw:53\n", mexprofrecord_[2]);
    } else if (strcmp(id, "*profile log*") == 0) {
        FILE* logfp;
        if (nrhs != 2 || mxGetString(prhs[1], id, sizeof(id)) != 0)
            mexErrMsgTxt("Must have two string arguments");
        logfp = fopen(id, "w+");
        if (!logfp)
            mexErrMsgTxt("Cannot open log for output");
        if (!mexprofrecord_)
            fprintf(logfp, "Profiler inactive\n");
        fprintf(logfp, "%d calls to gateway.mw:30\n", mexprofrecord_[1]);
        fprintf(logfp, "%d calls to gateway.mw:53\n", mexprofrecord_[2]);
        fclose(logfp);
    } else
        mexErrMsgTxt("Unknown identifier");
}
