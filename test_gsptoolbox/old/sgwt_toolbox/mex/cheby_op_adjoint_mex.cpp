#include "mex.h"
#include "cheby_op.h"

#define pma_y prhs[0]
#define pma_L prhs[1]
#define pma_c prhs[2]
#define pma_arange prhs[3]

char usage_string[]="call as cheby_op_adjoint_mex(y,L,c,arange)\n"
"or as cheby_op_adjoint_mex(y) to resuse previous L,c,arange\n";

/* persisent data for mex file : represented as global variables */
static int initialized=0;
static cheby_op *my_cheby_op;

void cleanup(){
  if(initialized){
    delete my_cheby_op;
    initialized=0;
  }
}

extern "C" void mexFunction(int nlhs,mxArray *plhs[],int nrhs,
                            const mxArray *prhs[]){
  if (!(nrhs==1 || nrhs==4)){
    mexPrintf("%s",usage_string);
    return;
  }
  //     Check validity of inputs.   
  if (!mxIsDouble(pma_y) || mxIsSparse(pma_y)){
    mexPrintf("%s",usage_string);
    return;
  }
  if (nrhs==4){
    if (!cheby_op::is_valid_input(pma_L,pma_c,pma_arange) ){
      mexPrintf("%s",usage_string);
      return;
    }
  }

  //  Handle persistent memory

  if (nrhs==4){
    cleanup();
    my_cheby_op=new cheby_op(pma_L,pma_c,pma_arange);
    mexMakeMemoryPersistent(my_cheby_op);
    mexAtExit(cleanup);
    initialized=1;
  }
  /* at this point, static data should be initialized. If not, then function was called first with only
     one argument, which is an error */
  if (!initialized)
    mexErrMsgTxt("not initialized : call with L,c,arange at least once first\n");
  
  int N=my_cheby_op->N;
  int s=my_cheby_op->s;
  
  if (mxGetNumberOfElements(pma_y)!=N*s)
    mexErrMsgTxt("mismatch between dimensions of y and size of SGWT transform\n");
  /* now, actually do computation */
  plhs[0]=mxCreateDoubleMatrix(N,1,mxREAL);

  // set up array of pointers for SGWT coefficients y
  double *y[s];
  for (int ks=0;ks<s;ks++)
    y[ks]=mxGetPr(pma_y)+ks*N;

  my_cheby_op->cheby_op_adjoint(mxGetPr(plhs[0]),y);
  //sgwt_cheby_op_adjoint(mxGetPr(plhs[0]),y,L,c,arange);
}
