#include "mex.h"
#include "cheby_op.h"

#define pma_f prhs[0]
#define pma_L prhs[1]
#define pma_c prhs[2]
#define pma_arange prhs[3]

char usage_string[]="call as cheby_op_mex(y,L,c,arange)\n"
"or as cheby_op_mex(f) to resuse previous L,c,arange\n";

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
  if (!mxIsDouble(pma_f) || mxIsSparse(pma_f)){
    mexPrintf("%s",usage_string);
    return;
  }
  if (nrhs==4){
    if (!cheby_op::is_valid_input(pma_L,pma_c,pma_arange) ){
      mexPrintf("%s",usage_string);
      return;
    }
  }

  // Handle persistent memory

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

  if (mxGetNumberOfElements(pma_f)!=N)
    mexErrMsgTxt("mismatch between dimensions of f and L\n");


  /* now, actually do computation */
  plhs[0]=mxCreateDoubleMatrix(N * s,1,mxREAL);
  // set up array of pointers, each of size N corresponding to
  // a column of res
  double *res[s];
  for (int ks=0;ks<s;ks++){
    res[ks]=mxGetPr(plhs[0])+ks*N;
  }
  my_cheby_op->cheby_op_forward(res,mxGetPr(pma_f));
  //sgwt_cheby_op(res,mxGetPr(pma_f),L,c,arange);
}
  

     
    
