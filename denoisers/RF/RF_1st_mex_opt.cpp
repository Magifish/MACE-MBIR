////////////////////////////////////////////////////////
// RF_1st_mex.cpp
// .mex version of 1st order recursive filter 
//
// xhat = RF_1st_mex(y, sigma_s, sigma_r, noise_sigma, sigma_g, max_itr, ref);
//
////////////////////////////////////////////////////////
#include <iostream> // for debugging
#include <fstream> // for debugging
#include <cstring> // array copy
#include <ctime>
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

void RF_1st_p1_main(double*, double*, double*, double, int, int);
void RF_1st_p2_main(double*, double*, double*, double*, double, double, int, int, int);

/////////////////////////
// RF Main Function
/////////////////////////
void RF_1st_p1_main(double *dIdx, double *dIdy, double *ref, 
        double sigma_g, int M, int N)
{
    int i, j, k;
    
    // Find dIdx and dIdy
    // Note: C++ image is represented as a flipped MATLAB image matrix
    dIdx[0] = 0; dIdy[0] = 0;
    for (j = 1; j < M; j++) {
        dIdy[j] = abs(ref[j] - ref[j-1]);
        dIdx[j] = 0;
    }
    for (i = M; i < (N-1)*M+1; i += M) {
        dIdx[i] = abs(ref[i] - ref[i-M]);
        dIdy[i] = 0;
    }
    for (i = 0; i < N-1; i++) {
        for (j = 0; j < M-1; j++) {
            k = (i+1)*M+(j+1);
            dIdy[k] = abs(ref[k] - ref[k-1]);
            dIdx[k] = abs(ref[k] - ref[k-M]);
        }
    }
}

void RF_1st_p2_main(double *xhat_, double *y_, double *dIdx_, double *dIdy_,
        double sigma_s, double sigma_r, double noise_sigma, int max_itr, int M, int N)
{
    int i, j, itr;
    double *xhat, *y, *dIdx, *dIdy;
    double *dHdx, *dVdy;
    double *alpha, *den, *num;
    double sigma_H_i, a;
//     std::ofstream myfile;
//     
//     myfile.open("log.txt", std::ios_base::app);
    
    // Flip matrices [0.005s]
    // ToDo: validate M and N for non-square inputs
    // ToDo: is flipping really necessary? remove flipping for shorter runtime
    xhat = new double[M*N];
    y    = new double[M*N];
    dIdx = new double[M*N];
    dIdy = new double[M*N];
    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            y[j*N+i]    = y_[i*M+j];
            dIdx[j*N+i] = dIdx_[i*M+j];
            dIdy[j*N+i] = dIdy_[i*M+j];
        }
    }
    //i = M; M = N; N = i; // also flip #row and #col
    
    // dHdx, dVdy [0.003s]
    dHdx = new double[M*N];
    dVdy = new double[M*N];
    for (i = 0; i < M; i++){
        for (j = 0; j < N; j++){
            xhat[i*N+j] = y[i*N+j];
            dHdx[i*N+j] = 1 + sigma_s/sigma_r * MAX(abs(dIdx[i*N+j])-noise_sigma,0);
            dVdy[i*N+j] = 1 + sigma_s/sigma_r * MAX(abs(dIdy[i*N+j])-noise_sigma,0);
        }
    }
    
    // Recursive filter iterations [0.025s]
    alpha = new double[M*N];
    num = new double[M*N];
    den = new double[M*N];
    for (itr = 0; itr < max_itr; itr++) {
        sigma_H_i = sigma_s * sqrt((double)3.0) * pow(2,(double)(max_itr-(itr+1))) / sqrt(pow(4,(double)max_itr)-1);
        a = exp(-sqrt((double)2.0) / sigma_H_i);
        
        // Horizontal
        //std::clock_t t1 = std::clock();
        for (i = 0; i < M; i++){
            for (j = 0; j < N; j++){
                num[j*M+i] = -exp(dHdx[i*N+j] * log(a));  // flipping
                den[j*M+i] = xhat[i*N+j];           // flipping
            }
        }
        i = M; M = N; N = i; // also flip #row and #col
        for (i = 0; i < M; i++){
            for (j = 0; j < N; j++){
                alpha[i*N+j] = num[i*N+j];
                num[i*N+j] = den[i*N+j];
                den[i*N+j] = 1.0;
            }
        }
        for (i = 1; i < M; i++) {
            for (j = 0; j < N; j++) {
                num[i*N+j] = num[i*N+j] - alpha[i*N+j]*num[(i-1)*N+j];
                den[i*N+j] = den[i*N+j] - alpha[i*N+j]*den[(i-1)*N+j];
            }
        }
        for (i = M-2; i >= 0; i--) {
            for (j = 0; j < N; j++) {
                num[i*N+j] = num[i*N+j] - alpha[(i+1)*N+j]*num[(i+1)*N+j];
                den[i*N+j] = den[i*N+j] - alpha[(i+1)*N+j]*den[(i+1)*N+j];
            }
        }
        for (i = 0; i < M; i++)
            for (j = 0; j < N; j++)
                xhat[j*M+i] = num[i*N+j] / den[i*N+j];
        i = M; M = N; N = i; // also flip #row and #col
        // Vertical
        for (i = 0; i < M; i++) {
            for (j = 0; j < N; j++) {
                alpha[i*N+j] = -exp(dVdy[i*N+j] * log(a));
                num[i*N+j] = xhat[i*N+j];
                den[i*N+j] = 1.0;
            }
        }
        for (i = 1; i < M; i++) {
            for (j = 0; j < N; j++) {
                num[i*N+j] = num[i*N+j] - alpha[i*N+j]*num[(i-1)*N+j];
                den[i*N+j] = den[i*N+j] - alpha[i*N+j]*den[(i-1)*N+j];
            }
        }
        for (i = M-2; i >= 0; i--) {
            for (j = 0; j < N; j++) {
                num[i*N+j] = num[i*N+j] - alpha[(i+1)*N+j]*num[(i+1)*N+j];
                den[i*N+j] = den[i*N+j] - alpha[(i+1)*N+j]*den[(i+1)*N+j];
            }
        }
        for (i = 0; i < M; i++)
            for (j = 0; j < N; j++)
                xhat[i*N+j] = num[i*N+j] / den[i*N+j];
        
    }
    
    // Flip xhat [0.001s]
    for (i = 0; i < M; i++)
        for (j = 0; j < N; j++)
            xhat_[j*M+i] = xhat[i*N+j];
    
//     myfile.close();
    delete den; delete num; delete alpha;
    delete dVdy; delete dHdx;
    delete dIdy; delete dIdx;
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
    double *ref;
    double sigma_s, sigma_r, noise_sigma, sigma_g;
    int max_itr;
    
    mxArray *dIdx, *dIdy;
    mxArray *lhs[1], *rhs[3];
    mxArray *mxstr;
    mxArray *gf;
    double *iptr, *dptr;
    
    int M, N, Mr, Nr;
    int r, win;
    int i;
    bool ref_flag = false;
    
    //std::ofstream myfile;
    //myfile.open("tmp.txt");
    
    /////////////////////
    // Inputs
    /////////////////////
    // Verifying
    if (nrhs < 3 || nrhs > 7)
        mexErrMsgTxt("Number of ARGS must be between 3 and 7.");
    
    // Input 0
    y = mxGetPr(prhs[0]);
    M = mxGetM(prhs[0]);
    N = mxGetN(prhs[0]);
    
    //myfile << M << std::endl;
    //myfile << N << std::endl;
    
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
    
    ////////////////////
    // Outputs
    ////////////////////
    dIdx = mxCreateDoubleMatrix(M, N, mxREAL);
    dIdy = mxCreateDoubleMatrix(M, N, mxREAL);
    RF_1st_p1_main(mxGetPr(dIdx), mxGetPr(dIdy), ref, sigma_g, M, N);

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
    rhs[0] = dIdx;
    rhs[1] = gf;
    rhs[2] = mxCreateString("replicate");
    mexCallMATLAB(1, lhs, 3, rhs, "imfilter");
    mxDestroyArray(dIdx);
    dIdx = lhs[0];
    //dIdy = imfilter(dIdy, gf, 'replicate');
    rhs[0] = dIdy;
    mexCallMATLAB(1, lhs, 3, rhs, "imfilter");
    mxDestroyArray(dIdy);
    dIdy = lhs[0];
    mxDestroyArray(rhs[2]);
    mxDestroyArray(gf);
    //mexCallMATLAB(0, NULL, 1, &dIdy, "disp");
    
    plhs[0] = mxCreateDoubleMatrix(M, N, mxREAL);
    xhat    = mxGetPr(plhs[0]);
    
    RF_1st_p2_main(xhat, y, mxGetPr(dIdx), mxGetPr(dIdy), sigma_s, sigma_r, noise_sigma, max_itr, M, N);

    mxDestroyArray(dIdy);
    mxDestroyArray(dIdx);
    
    //myfile << M << std::endl;
    //myfile << N << std::endl;
    
    //myfile.close();

    
    if (ref_flag)
        delete ref;
}