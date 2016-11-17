#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>  //  for toupper() macro in streq()

#ifdef WIN32
#include <windows.h>
#endif

#include <cmath>    //  Definitions for Matlab API
#include <mex.h>    //  Definitions for Matlab API

using namespace std;


void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  double *in,*com,*Mnew;
  int debug = 0;

  /* Check for proper number of arguments. */
  if(nrhs<2) {
    mexErrMsgTxt("At least two input required.");
  } else if(nrhs>2) {
    mexErrMsgTxt("Too many output arguments");
  }
  
  /* The input must be a noncomplex scalar double.*/
  int N1 = mxGetM(prhs[0]);
  int M1 = mxGetN(prhs[0]);
  if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
      !(M1==N1) ) {
    mexErrMsgTxt("Input matrix must be square and non-complex.");
  }

  int N2 = mxGetM(prhs[1]);
  int M2 = mxGetN(prhs[1]);
  if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ||
      !(N2 == 1) || (M2 != M1)) {
    mexErrMsgTxt("Community should be a non-complex vector.");
  }
  
  int N = M2;

  /* Assign pointers to each input and output. */
  in = mxGetPr(prhs[0]);
  com = mxGetPr(prhs[1]);

  int M = 0;
  for (int i = 0; i < N; i++){
    int t = (int)floor(com[i]);
    if (t > M){
      M = t;
    }
  }
  
  /* Create matrix for the return argument. */
  plhs[0] = mxCreateDoubleMatrix(M,M,mxREAL);
  Mnew = mxGetPr(plhs[0]);

  for (int i = 0; i < M; i++){
    for (int j = 0; j < M; j++){
      Mnew[i*M+j] = 0.0;
    }
  }

  for (int i = 0; i < N; i++){
    for (int j = 0; j < N; j++){
      int comi = (int)floor(com[i])-1;
      int comj = (int)floor(com[j])-1;
      if (comi >= M || comi < 0 || comj >= M || comj < 0){
	mexErrMsgTxt("Community number error.");
      }
      Mnew[comi*M+comj] += in[i*N+j];
    }
  }
  
  
}

