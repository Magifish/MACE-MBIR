////////////////////////////////////////////////////////
// RF_3rd_mex.cpp
// .mex version of 3rd order recursive filter 
//
// xhat = RF_3rd_mex(y, sigma_s, sigma_r, noise_sigma, sigma_g, max_itr, ref, a);
//
////////////////////////////////////////////////////////
#include <math.h>
#include "mex.h"
#include "matrix.h"

#ifndef PI
#define PI 3.1415926
#endif
#ifndef MIN
#define	MIN(a, b) ((a) < (b) ? (a) : (b))
#endif
#ifndef MAX
#define	MAX(a, b) ((a) > (b) ? (a) : (b))
#endif

void RF_3rd_p1_main(double*, double*, double*, double*, double*, double*,
        double*, double, int, int);
void RF_3rd_p2_main(double*, double*, double*, double*, double*, double*,
        double*, double*, double, double, double, int, double*, int, int);

/////////////////////////
// RF Main Function
/////////////////////////
void RF_3rd_p1_main(double *dIdx1, double *dIdy1, double *dIdx2, double *dIdy2,
        double *dIdx3, double *dIdy3, double *ref, double sigma_g, int M, int N)
{
    int i, j, k;
    //double *ds[6] = {dIdx1, dIdy1, dIdx2, dIdy2, dIdx3, dIdy3};
    
    // Find dIdx's
    for (i = 0; i < N; i++) {
        dIdx1[i] = 0; dIdx2[i] = 0; dIdx3[i] = 0;
        dIdx1[i+M] = abs(ref[i+M] - ref[i]);
        dIdx2[i+M] = abs(ref[i+2*M] - ref[i]);
        dIdx3[i+M] = 0;
        k = i+(N-1)*M;
        dIdx1[k] = abs(ref[k] - ref[k-M]);
        dIdx2[k] = 0; dIdx3[k] = 0;
    }
    for (i = 2; i < N-1; i++) {
        for (j = 0; j < M; j++) {
            k = i*M+j;
            dIdx1[k] = abs(ref[k] - ref[k-M]);
            dIdx2[k] = abs(ref[k+M] - ref[k-M]);
            dIdx3[k] = abs(ref[k+M] - ref[k-2*M]);
        }
    }
    // Find dIdy's
    for (i = 0; i < M*N; i += M) {
        dIdy1[i] = 0; dIdy2[i] = 0; dIdy3[i] = 0;
        dIdy1[i+1] = abs(ref[i+1] - ref[i]);
        dIdy2[i+1] = abs(ref[i+2] - ref[i]);
        dIdy3[i+1] = 0;
        k = i+(M-1);
        dIdy1[k] = abs(ref[k] - ref[k-1]);
        dIdy2[i+2] = 0; dIdy3[i+2] = 0;
    }
    for (i = 0; i < N; i++) {
        for (j = 2; j < M-1; j++) {
            k = i*M+j;
            dIdy1[k] = abs(ref[k] - ref[k-1]);
            dIdy2[k] = abs(ref[k+1] - ref[k-1]);
            dIdy3[k] = abs(ref[k+1] - ref[k-2]);
        }
    }
    
}

void RF_3rd_p2_main(double *xhat_, double *y_, double *dIdx1_, double *dIdy1_,
        double *dIdx2_, double *dIdy2_, double *dIdx3_, double *dIdy3_, double sigma_s,
        double sigma_r, double noise_sigma, int max_itr, double *a, int M, int N)
{
    int i, j, itr;
    double *xhat, *y;
    double *dIdx1, *dIdy1, *dIdx2, *dIdy2, *dIdx3, *dIdy3;
    double *dHdx1, *dVdy1, *dHdx2, *dVdy2, *dHdx3, *dVdy3;
    double *alpha1, *alpha2, *alpha3, *den, *num;
    double *temp1, *temp2, *temp3;
    double sigma_H_i, a_param;
    
    // Flip matrices
    // ToDo: validate M and N for non-square inputs
    // ToDo: is flipping really necessary? remove flipping for shorter runtime
    xhat = new double[M*N];
    y    = new double[M*N];
    dIdx1 = new double[M*N]; dIdy1 = new double[M*N];
    dIdx2 = new double[M*N]; dIdy2 = new double[M*N];
    dIdx3 = new double[M*N]; dIdy3 = new double[M*N];
    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            y[j*N+i]    = y_[i*M+j];
            dIdx1[j*N+i] = dIdx1_[i*M+j];
            dIdy1[j*N+i] = dIdy1_[i*M+j];
            dIdx2[j*N+i] = dIdx2_[i*M+j];
            dIdy2[j*N+i] = dIdy2_[i*M+j];
            dIdx3[j*N+i] = dIdx3_[i*M+j];
            dIdy3[j*N+i] = dIdy3_[i*M+j];
        }
    }
    //i = M; M = N; N = i; // also flip #row and #col
    
    // dHdx, dVdy
    dHdx1 = new double[M*N]; dVdy1 = new double[M*N];
    dHdx2 = new double[M*N]; dVdy2 = new double[M*N];
    dHdx3 = new double[M*N]; dVdy3 = new double[M*N];
    for (i = 0; i < M; i++){
        for (j = 0; j < N; j++){
            xhat[i*N+j] = y[i*N+j];
            dHdx1[i*N+j] = 1 + sigma_s/sigma_r * MAX(abs(dIdx1[i*N+j])-noise_sigma,0);
            dVdy1[i*N+j] = 1 + sigma_s/sigma_r * MAX(abs(dIdy1[i*N+j])-noise_sigma,0);
            dHdx2[i*N+j] = 2 + sigma_s/sigma_r * MAX(abs(dIdx2[i*N+j])-noise_sigma,0);
            dVdy2[i*N+j] = 2 + sigma_s/sigma_r * MAX(abs(dIdy2[i*N+j])-noise_sigma,0);
            dHdx3[i*N+j] = 3 + sigma_s/sigma_r * MAX(abs(dIdx3[i*N+j])-noise_sigma,0);
            dVdy3[i*N+j] = 3 + sigma_s/sigma_r * MAX(abs(dIdy3[i*N+j])-noise_sigma,0);
        }
    }

    // Recursive filter iterations
    alpha1 = new double[M*N];
    alpha2 = new double[M*N];
    alpha3 = new double[M*N];
    temp1 = new double[M*N];
    temp2 = new double[M*N];
    temp3 = new double[M*N];
    num = new double[M*N];
    den = new double[M*N];
    for (itr = 0; itr < max_itr; itr++) {
        sigma_H_i = sigma_s * sqrt((double)3.0) * pow(2,(double)(max_itr-(itr+1))) / sqrt(pow(4,(double)max_itr)-1);
        a_param = exp(-sqrt((double)2.0) / sigma_H_i);
        // Horizontal
        for (i = 0; i < M; i++) {
            for (j = 0; j < N; j++) {
                temp1[j*M+i] = exp(dHdx1[i*N+j] * log(a_param));  // flipping
                temp2[j*M+i] = exp(dHdx2[i*N+j] * log(a_param));  // flipping
                temp3[j*M+i] = exp(dHdx3[i*N+j] * log(a_param));  // flipping
                den[j*M+i] = xhat[i*N+j];                   // flipping
            }
        }
        i = M; M = N; N = i; // also flip #row and #col
        for (i = 0; i < M; i++) {
            for (j = 0; j < N; j++) {
                alpha1[i*N+j] = temp1[i*N+j];
                alpha2[i*N+j] = temp2[i*N+j];
                alpha3[i*N+j] = temp3[i*N+j];
                num[i*N+j] = den[i*N+j];
                den[i*N+j] = 1.0;
            }
        }
        for (i = 1; i < M; i++) {
            for (j = 0; j < N; j++) {
                //num[i*N+j] = num[i*N+j] - a[0]*alpha1[i*N+j]*num[(i-1)*N+j];
                //den[i*N+j] = den[i*N+j] - a[0]*alpha1[i*N+j]*den[(i-1)*N+j];
                if (i == 1) {
                    num[i*N+j] = num[i*N+j] - a[0]*alpha1[i*N+j]*num[(i-1)*N+j];
                    den[i*N+j] = den[i*N+j] - a[0]*alpha1[i*N+j]*den[(i-1)*N+j];
                } else if (i == 2) {
                    num[i*N+j] = num[i*N+j] - a[0]*alpha1[i*N+j]*num[(i-1)*N+j] - a[1]*alpha2[i*N+j]*num[(i-2)*N+j];
                    den[i*N+j] = den[i*N+j] - a[0]*alpha1[i*N+j]*den[(i-1)*N+j] - a[1]*alpha2[i*N+j]*den[(i-2)*N+j];
                } else {
                    num[i*N+j] = num[i*N+j] - a[0]*alpha1[i*N+j]*num[(i-1)*N+j] - a[1]*alpha2[i*N+j]*num[(i-2)*N+j] - a[2]*alpha3[i*N+j]*num[(i-3)*N+j];
                    den[i*N+j] = den[i*N+j] - a[0]*alpha1[i*N+j]*den[(i-1)*N+j] - a[1]*alpha2[i*N+j]*den[(i-2)*N+j] - a[2]*alpha3[i*N+j]*den[(i-3)*N+j];
                }
            }
        }
        for (i = M-2; i >= 0; i--) {
            for (j = 0; j < N; j++) {
                //num[i*N+j] = num[i*N+j] - a[0]*alpha1[(i+1)*N+j]*num[(i+1)*N+j];
                //den[i*N+j] = den[i*N+j] - a[0]*alpha1[(i+1)*N+j]*den[(i+1)*N+j];
                if (i == M-2) {
                    num[i*N+j] = num[i*N+j] - a[0]*alpha1[(i+1)*N+j]*num[(i+1)*N+j];
                    den[i*N+j] = den[i*N+j] - a[0]*alpha1[(i+1)*N+j]*den[(i+1)*N+j];
                } else if (i == M-3) {
                    num[i*N+j] = num[i*N+j] - a[0]*alpha1[(i+1)*N+j]*num[(i+1)*N+j] - a[1]*alpha2[(i+1)*N+j]*num[(i+2)*N+j];
                    den[i*N+j] = den[i*N+j] - a[0]*alpha1[(i+1)*N+j]*den[(i+1)*N+j] - a[1]*alpha2[(i+1)*N+j]*den[(i+2)*N+j];
                } else {
                    num[i*N+j] = num[i*N+j] - a[0]*alpha1[(i+1)*N+j]*num[(i+1)*N+j] - a[1]*alpha2[(i+1)*N+j]*num[(i+2)*N+j] - a[2]*alpha3[(i+1)*N+j]*num[(i+3)*N+j];
                    den[i*N+j] = den[i*N+j] - a[0]*alpha1[(i+1)*N+j]*den[(i+1)*N+j] - a[1]*alpha2[(i+1)*N+j]*den[(i+2)*N+j] - a[2]*alpha3[(i+1)*N+j]*den[(i+3)*N+j];
                }
            }
        }
        for (i = 0; i < M; i++)
            for (j = 0; j < N; j++)
                xhat[j*M+i] = num[i*N+j] / den[i*N+j];
        i = M; M = N; N = i; // also flip #row and #col
        // Vertical
        for (i = 0; i < M; i++) {
            for (j = 0; j < N; j++) {
                alpha1[i*N+j] = exp(dVdy1[i*N+j] * log(a_param));
                alpha2[i*N+j] = exp(dVdy2[i*N+j] * log(a_param));
                alpha3[i*N+j] = exp(dVdy3[i*N+j] * log(a_param));
                num[i*N+j] = xhat[i*N+j];
                den[i*N+j] = 1.0;
            }
        }
        for (i = 1; i < M; i++) {
            for (j = 0; j < N; j++) {
//                 num[i*N+j] = num[i*N+j] - a[0]*alpha1[i*N+j]*num[(i-1)*N+j];
//                 den[i*N+j] = den[i*N+j] - a[0]*alpha1[i*N+j]*den[(i-1)*N+j];
                if (i == 1) {
                    num[i*N+j] = num[i*N+j] - a[0]*alpha1[i*N+j]*num[(i-1)*N+j];
                    den[i*N+j] = den[i*N+j] - a[0]*alpha1[i*N+j]*den[(i-1)*N+j];
                } else if (i == 2) {
                    num[i*N+j] = num[i*N+j] - a[0]*alpha1[i*N+j]*num[(i-1)*N+j] - a[1]*alpha2[i*N+j]*num[(i-2)*N+j];
                    den[i*N+j] = den[i*N+j] - a[0]*alpha1[i*N+j]*den[(i-1)*N+j] - a[1]*alpha2[i*N+j]*den[(i-2)*N+j];
                } else {
                    num[i*N+j] = num[i*N+j] - a[0]*alpha1[i*N+j]*num[(i-1)*N+j] - a[1]*alpha2[i*N+j]*num[(i-2)*N+j] - a[2]*alpha3[i*N+j]*num[(i-3)*N+j];
                    den[i*N+j] = den[i*N+j] - a[0]*alpha1[i*N+j]*den[(i-1)*N+j] - a[1]*alpha2[i*N+j]*den[(i-2)*N+j] - a[2]*alpha3[i*N+j]*den[(i-3)*N+j];
                }
            }
        }
        for (i = M-2; i >= 0; i--) {
            for (j = 0; j < N; j++) {
//                 num[i*N+j] = num[i*N+j] - a[0]*alpha1[(i+1)*N+j]*num[(i+1)*N+j];
//                 den[i*N+j] = den[i*N+j] - a[0]*alpha1[(i+1)*N+j]*den[(i+1)*N+j];
                if (i == M-2) {
                    num[i*N+j] = num[i*N+j] - a[0]*alpha1[(i+1)*N+j]*num[(i+1)*N+j];
                    den[i*N+j] = den[i*N+j] - a[0]*alpha1[(i+1)*N+j]*den[(i+1)*N+j];
                } else if (i == M-3) {
                    num[i*N+j] = num[i*N+j] - a[0]*alpha1[(i+1)*N+j]*num[(i+1)*N+j] - a[1]*alpha2[(i+1)*N+j]*num[(i+2)*N+j];
                    den[i*N+j] = den[i*N+j] - a[0]*alpha1[(i+1)*N+j]*den[(i+1)*N+j] - a[1]*alpha2[(i+1)*N+j]*den[(i+2)*N+j];
                } else {
                    num[i*N+j] = num[i*N+j] - a[0]*alpha1[(i+1)*N+j]*num[(i+1)*N+j] - a[1]*alpha2[(i+1)*N+j]*num[(i+2)*N+j] - a[2]*alpha3[(i+1)*N+j]*num[(i+3)*N+j];
                    den[i*N+j] = den[i*N+j] - a[0]*alpha1[(i+1)*N+j]*den[(i+1)*N+j] - a[1]*alpha2[(i+1)*N+j]*den[(i+2)*N+j] - a[2]*alpha3[(i+1)*N+j]*den[(i+3)*N+j];
                }
            }
        }
        for (i = 0; i < M; i++)
            for (j = 0; j < N; j++)
                xhat[i*N+j] = num[i*N+j] / den[i*N+j];
    }
    
    // Flip xhat
    for (i = 0; i < M; i++)
        for (j = 0; j < N; j++)
            xhat_[j*M+i] = xhat[i*N+j];
    
    delete den; delete num;
    delete temp3; delete temp2; delete temp1;
    delete alpha3; delete alpha2; delete alpha1;
    delete dVdy3; delete dHdx3;
    delete dVdy2; delete dHdx2;
    delete dVdy1; delete dHdx1;
    delete dIdy3; delete dIdx3;
    delete dIdy2; delete dIdx2;
    delete dIdy1; delete dIdx1;
    delete y; delete xhat;
}


/////////////////////////
// Interface Function
/////////////////////////
void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    // Variable Declarations
    double *xhat;
    double *y;
    double *ref, *a;
    double sigma_s, sigma_r, noise_sigma, sigma_g;
    int max_itr;
    
    mxArray *dIdx1, *dIdy1, *dIdx2, *dIdy2, *dIdx3, *dIdy3;
    mxArray *lhs[1], *rhs[3];
    mxArray *mxstr;
    mxArray *gf;
    double *iptr, *dptr;
    
    int M, N, Mr, Nr;
    int r, win;
    int i;
    int a_len;
    bool ref_flag = false;
    
    /////////////////////
    // Inputs
    /////////////////////
    // Verifying
    if (nrhs < 3 || nrhs > 8)
        mexErrMsgTxt("Number of ARGS must be between 3 and 8.");
    
    // Input 0
    y = mxGetPr(prhs[0]);
    M = mxGetM(prhs[0]);
    N = mxGetN(prhs[0]);
    
    // Input 1
    sigma_s = mxGetScalar(prhs[1]);
    
    // Input 2
    sigma_r = mxGetScalar(prhs[2]);
    
    // Input 3
    if (nrhs < 4) noise_sigma = 0;
    else noise_sigma = mxGetScalar(prhs[3]);
    
    // Input 4
    if (nrhs < 5) sigma_g = 0.5;
    else sigma_g = mxGetScalar(prhs[4]);
    
    // Input 5
    if (nrhs < 6) max_itr = 3;
    else max_itr = mxGetScalar(prhs[5]);
    
    // Input 6
    if (nrhs < 7) {
        ref = new double[M*N];
        ref_flag = true;
        for (i = 0; i < M*N; i++) {
            ref[i] = y[i];
        }
    } else {
        ref = mxGetPr(prhs[6]);
        Mr = mxGetM(prhs[6]);
        Nr = mxGetN(prhs[6]);
        if (Mr != M || Nr != N)
            mexErrMsgTxt("Noisy image size must match reference image size.");
    }
    
    // Input 7
    if (nrhs < 8) {
        a = new double[3];
        a[0] = -1; a[1] = 0; a[2] = 0;
        a_len = -1;
    } else {
        a = mxGetPr(prhs[7]);
        a_len = mxGetN(prhs[7]);
        if (a_len != 3)
            mexErrMsgTxt("a array length must be 3.");
    }
    
    ////////////////////
    // Outputs
    ////////////////////
    dIdx1 = mxCreateDoubleMatrix(M, N, mxREAL);
    dIdy1 = mxCreateDoubleMatrix(M, N, mxREAL);
    dIdx2 = mxCreateDoubleMatrix(M, N, mxREAL);
    dIdy2 = mxCreateDoubleMatrix(M, N, mxREAL);
    dIdx3 = mxCreateDoubleMatrix(M, N, mxREAL);
    dIdy3 = mxCreateDoubleMatrix(M, N, mxREAL);
    RF_3rd_p1_main(mxGetPr(dIdx1), mxGetPr(dIdy1), mxGetPr(dIdx2), mxGetPr(dIdy2),
            mxGetPr(dIdx3), mxGetPr(dIdy3), ref, sigma_g, M, N);
    
    r = (int)MAX(ceil(sqrt(-2*pow(sigma_g,2)*log(0.001*sqrt(2*PI*pow(sigma_g,2))))),1);
    win = r*2+1;
    
    //gf = fspecial('gaussian',[win win],sigma_g);
    rhs[0] = mxCreateString("gaussian");
    rhs[1] = mxCreateNumericMatrix(1, 2, mxDOUBLE_CLASS, mxREAL);
    rhs[2] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
    iptr = mxGetPr(rhs[1]);
    iptr[0] = win; iptr[1] = win;
    dptr = mxGetPr(rhs[2]);
    dptr[0] = sigma_g;
    mexCallMATLAB(1, lhs, 3, rhs, "fspecial");
    mxDestroyArray(rhs[2]);
    mxDestroyArray(rhs[1]);
    mxDestroyArray(rhs[0]);
    gf = lhs[0];
    //dIdx = imfilter(dIdx, gf, 'replicate');
    rhs[0] = dIdx1;
    rhs[1] = gf;
    rhs[2] = mxCreateString("replicate");
    mexCallMATLAB(1, lhs, 3, rhs, "imfilter");
    mxDestroyArray(dIdx1);
    dIdx1 = lhs[0];
    //dIdy = imfilter(dIdy, gf, 'replicate');
    rhs[0] = dIdy1;
    mexCallMATLAB(1, lhs, 3, rhs, "imfilter");
    mxDestroyArray(dIdy1);
    dIdy1 = lhs[0];
    //dIdy = imfilter(dIdx, gf, 'replicate');
    rhs[0] = dIdx2;
    mexCallMATLAB(1, lhs, 3, rhs, "imfilter");
    mxDestroyArray(dIdx2);
    dIdx2 = lhs[0];
    //dIdy = imfilter(dIdy, gf, 'replicate');
    rhs[0] = dIdy2;
    mexCallMATLAB(1, lhs, 3, rhs, "imfilter");
    mxDestroyArray(dIdy2);
    dIdy2 = lhs[0];
    //dIdy = imfilter(dIdx, gf, 'replicate');
    rhs[0] = dIdx3;
    mexCallMATLAB(1, lhs, 3, rhs, "imfilter");
    mxDestroyArray(dIdx3);
    dIdx3 = lhs[0];
    //dIdy = imfilter(dIdy, gf, 'replicate');
    rhs[0] = dIdy3;
    mexCallMATLAB(1, lhs, 3, rhs, "imfilter");
    mxDestroyArray(dIdy3);
    dIdy3 = lhs[0];
    mxDestroyArray(rhs[2]);
    mxDestroyArray(gf);
    //mexCallMATLAB(0, NULL, 1, &dIdy, "disp");
    
    plhs[0] = mxCreateDoubleMatrix(M, N, mxREAL);
    xhat    = mxGetPr(plhs[0]);
    
    RF_3rd_p2_main(xhat, y, mxGetPr(dIdx1), mxGetPr(dIdy1), mxGetPr(dIdx2),
            mxGetPr(dIdy2), mxGetPr(dIdx3), mxGetPr(dIdy3), 
            sigma_s, sigma_r, noise_sigma, max_itr, a, M, N);
    
    mxDestroyArray(dIdy3);
    mxDestroyArray(dIdx3);
    mxDestroyArray(dIdy2);
    mxDestroyArray(dIdx2);
    mxDestroyArray(dIdy1);
    mxDestroyArray(dIdx1);
    
    if (a_len == -1) delete a;
    if (ref_flag) delete ref;
}