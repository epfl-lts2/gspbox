// Implementation : Antoine Scherrer
// antoine.scherrer@ens-lyon.fr
// After "" 
// "Fast unfolding of community hierarchies in large networks"
// Vincent D. Blondel, Jean-Loup Guillaume, Renaud Lambiotte,
// Etienne Lefebvre
// http://arxiv.org/abs/0803.0476

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

int find_comty(double out[], double in[],int N,int debug)
{

  vector<double>DeltaQ;
  vector<int>COM;
  vector<int>CSize;
  vector<double>K;
  vector<double>Cost;
  vector<double>SumIn;
  vector<double>SumTot;
  vector<double>Self;

  double min_increase = 0.0;
  int max_pass = 100;

  for (int i = 0;i < N; i++){
    SumIn.push_back(in[i*N+i]);
    SumTot.push_back(0.0);
    DeltaQ.push_back(0.0);
    K.push_back(0.0);
    Cost.push_back(0.0);
    COM.push_back(i);
    Self.push_back(0.0);
    CSize.push_back(1);
  }

  double m = 0.0;
  double m2 = 0.0;
  for (int i = 0;i < N; i++){
    double mt = 0.0;
    for (int j = 0;j < N; j++){
      m2 += in[i*N+j];
      mt += in[i*N+j];
      if (i==j)
	Self[i] = in[i*N+j];
    }
    K[i] = mt;
    SumTot[i] = mt;
  }
  m = m2 / 2;

  bool gain = true;
  int Niter = 0;

  double sum_g = 0.0;
  double mod_exp = 0.0;

  double init_mod = 0.0;
  for (int k = 0; k < N; k++){
    init_mod += SumIn[k]/m2 - (SumTot[k]/m2)*(SumTot[k]/m2);
  }
  if (debug){
    mexPrintf("Inital Mod=%e \n",init_mod);
  }
  double new_mod = init_mod;
  double cur_mod = init_mod;

  // Main loop 
  while (gain && Niter < max_pass){

    cur_mod = new_mod;
    sum_g = 0.0;
    mod_exp = 0.0;
    gain = false;

    // Loop over all nodes
    for (int i = 0; i < N; i++){
      int Ci = COM[i];
      double best_increase = 0.0;
      int best_com = Ci;
      for (int l = 0; l < N; l++){
	DeltaQ[l] = 0.0;
      }
            
      //Delete i from its community
      COM[i] = -1;
      for (int k = 0; k < N; k++){
	if (COM[k] == Ci){	   
	  SumIn[Ci] -= 2.0*in[i*N+k];
	}
      }
      SumIn[Ci] -= Self[i];
      SumTot[Ci] -= K[i];
      CSize[Ci]--;

      //Loop over neightbor
      for (int j = 0; j < N; j++){
	//Check if neighbor and different com
	int Cj = COM[j];
	if (in[i*N+j] != 0 && DeltaQ[Cj] == 0 && Cj != -1){
	  //Compute Ki_in
	  double Ki_in = 0.0;
	  for (int k = 0; k < N; k++){
	    if (COM[k] == Cj){
	      Ki_in += 2.0*in[i*N+k];
	    }
	  }
	  //DeltaQ[Cj] = (SumIn[Cj]+2.0*Ki_in)/m2 - ((SumTot[Cj]+K[i])/m2)*((SumTot[Cj]+K[i])/m2) - ( SumIn[Cj]/m2 - (SumTot[Cj]/m2)*(SumTot[Cj]/m2) - (K[i]/m2)*(K[i]/m2));
	  DeltaQ[Cj] = (Ki_in)/m2 - (2.0 * K[i] * SumTot[Cj]) / (m2*m2);
	  if (debug){
	    //mexPrintf("Gain for comm %d => %f\n",Cj,DeltaQ[Cj]);
	  }
	  if (DeltaQ[Cj] > best_increase){
	    best_increase = DeltaQ[Cj];
	    best_com = Cj;
	  }
	}
      }
      
      if (best_increase > 0){
	// Move i in highest gain community
	if (best_com < 0 || best_com > N-1){
	  mexPrintf("Internal ERROR (%d)\n",best_com);
	}
	if (debug){
	  mexPrintf("Move %d => %d\n",i,best_com);
	}
	sum_g += best_increase;
	Cost[i] = best_increase;
	mod_exp += best_increase;

	if (best_com != Ci){
	  gain = true;
	}
      }  else {
	best_com = Ci;
      }

      // Update com
      for (int k = 0; k < N; k++){
	if (COM[k] == best_com){
	  SumIn[best_com] += 2.0*in[i*N+k];
	} 
      }
      SumIn[best_com] += Self[i];
      SumTot[best_com] += K[i];
      COM[i] = best_com;
      CSize[best_com]++;
                  
    }
    
    new_mod = 0.0;
    for (int k = 0; k < N; k++){
      if (SumTot[k] > 0) {
	new_mod += SumIn[k]/m2 - (SumTot[k]/m2)*(SumTot[k]/m2);
      }
    }

    if (new_mod-cur_mod <= min_increase){
      gain = false;
    }
    
    if (debug){
      int ncom1 = 0;
      int ncom2 = 0;
      for (int k = 0; k < N; k++){
	if (CSize[k] > 0){
	  ncom1++;
	}
	if (CSize[k] > 1){
	  ncom2++;
	}
      }
      mexPrintf("It %d - Mod=%f (%f - %e) %d com (%d non isolated)\n",Niter,new_mod,cur_mod,new_mod-cur_mod,ncom1,ncom2);
    }

    Niter++;

  }

  //return result
  for (int k = 0; k < N; k++){
    out[k] = COM[k]+1;
  }
  
  return Niter;
}

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  double *in,*com,*deb,*niter;
  int debug = 0;

  /* Check for proper number of arguments. */
  if(nrhs<1) {
    mexErrMsgTxt("At least One input required.");
  } else if(nrhs==1) {
    debug = 0;
  } else if(nrhs==2) {
    deb = mxGetPr(prhs[1]);
    debug = (int)ceil(deb[0]);
  } else if(nrhs>2) {
    mexErrMsgTxt("Too many input arguments");
  }

  /* Check for proper number of arguments. */
  if(nlhs<2) {
    mexErrMsgTxt("Two output required.");
  } else if(nrhs>2) {
    mexErrMsgTxt("Too many output arguments");
  }
  
  /* The input must be a noncomplex scalar double.*/
  int N = mxGetM(prhs[0]);
  int M = mxGetN(prhs[0]);
  if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
      !(M==N) ) {
    mexErrMsgTxt("Input matrix must be square and non-complex.");
  }
  
  /* Create matrix for the return argument. */
  plhs[0] = mxCreateDoubleMatrix(1,N, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(1,1, mxREAL);
  
  /* Assign pointers to each input and output. */
  in = mxGetPr(prhs[0]);
  com = mxGetPr(plhs[0]);
  niter = mxGetPr(plhs[1]);
  
  /* Call the timestwo subroutine. */
  //find_comty(out,in,N);
  niter[0] = (double)find_comty(com,in,N,debug);
  
}

