#include <exception>
#include <algorithm>
#include <stdexcept>
#include "mex.h"
#include "cheby_op.h"

using std::copy;

void sparse_compressed_row::sparse_matrix_multiply(const double *prx,double *res){
  double tmp;
  for (int m=0;m<M;m++){
    tmp=0;
    for(mwIndex k=0;k<jc[m+1]-jc[m];k++)
      tmp+=prx[ ir[jc[m]+k] ] * pr[ jc[m]+k ]; 
    res[m]=tmp;
  }
}

sparse_compressed_row::sparse_compressed_row(const mxArray *array_in){
  if (!is_valid_input(array_in))
    throw std::runtime_error("matrix must be sparse and of double type");
  M=mxGetM(array_in);
  N=mxGetN(array_in);
  // We want to store compressed-row form of matrix. MATLAB uses compressed
  // column form. The compressed-column form for the transpose can be interpreted
  // as the compressed-row form of the original matrix. 
  mxArray *arg[1];
  mxArray *transposed[1];
  arg[0]=mxDuplicateArray(array_in);
  mexCallMATLAB(1,transposed,1,arg,"transpose");
  // copy memory, as returned array will not be persistent
  nz=mxGetNzmax(transposed[0]);
  pr=new double[nz];
  ir=new mwIndex[nz];
  jc=new mwIndex[M+1];
  copy(mxGetPr(transposed[0]),mxGetPr(transposed[0])+nz  ,pr);
  copy(mxGetIr(transposed[0]),mxGetIr(transposed[0])+nz  ,ir);
  copy(mxGetJc(transposed[0]),mxGetJc(transposed[0])+M+1 ,jc);  
}
sparse_compressed_row::~sparse_compressed_row(){
  delete pr;
  delete ir;
  delete jc;
}
bool sparse_compressed_row::is_valid_input(const mxArray *array_in){
  return mxIsSparse(array_in) && mxIsDouble(array_in);
}


cheby_coeff::cheby_coeff(const mxArray *c_in){
  if (!is_valid_input(c_in)){
    throw std::runtime_error("cheby_coeff constructor called with invalid input");
  }
  s=mxGetNumberOfElements(c_in);
  m=mxGetNumberOfElements(mxGetCell(c_in,0));
  c.resize(s);
  double *pr_tmp;
  for (int ks=0;ks<s;ks++){
    c[ks].resize(m);
    pr_tmp=mxGetPr(mxGetCell(c_in,ks));
    copy(pr_tmp,pr_tmp+m,c[ks].begin());
  }
}
// this constructor copies a single polynomial from another cheby_coeff
// this is useful for adjoint computation
// ks must be in range 0<=ks < c_in->s
cheby_coeff::cheby_coeff(cheby_coeff *c_in,int ks){
  s=1;
  m=c_in->m;
  c.resize(1);
  c[0].resize(m);
  copy(c_in->c[ks].begin(),c_in->c[ks].end(),c[0].begin());
}

//cheby_coeff::~cheby_coeff(){}

bool cheby_coeff::is_valid_input(const mxArray *c_in){
  // check that array is cell, and that every cell is double type
  // with same # of elements
  if (!mxIsCell(c_in))
    return false;
  unsigned int s=mxGetNumberOfElements(c_in);
  unsigned int m=mxGetNumberOfElements(mxGetCell(c_in,0));
  bool is_valid=true;
  mxArray *tmp;
  for (unsigned int k=0;k<s;k++){
    tmp=mxGetCell(c_in,k);
    if (!mxIsDouble(tmp) || mxGetNumberOfElements(tmp)!=m){
      is_valid=false;
      break;
    }
  }
  return is_valid;
}


approximation_range::approximation_range(const mxArray *array_in){
  if (!is_valid_input(array_in) ){
    throw std::runtime_error("invalid input for approximation_range constructor");
  }
  arange_[0]=mxGetPr(array_in)[0];
  arange_[1]=mxGetPr(array_in)[1];  
}

double approximation_range::operator[](int k){
  return arange_[k];
}

bool approximation_range::is_valid_input(const mxArray *array_in){
  return mxIsDouble(array_in) && mxGetNumberOfElements(array_in)==2 && !mxIsSparse(array_in);
}

// forward computation
// actually do the computation. res[k] has L->N elements, f has L->N elements
void sgwt_cheby_op(double *res[],double *f,sparse_compressed_row *L, cheby_coeff *c, approximation_range *arange){
  double a1=((*arange)[1]-(*arange)[0])/2;
  double a2=((*arange)[1]+(*arange)[0])/2;
  int N=L->N;
  // later, these should be made class members, so that allocation/deallocation
  // occurs whenever static L,c,arange are updated. For now, just allocate/ deallocate inside this 
  // function
  double *Twf_old ;
  double *Twf_cur ;
  double *Twf_new ;
  double *tmp     ;
  Twf_old = new double[N];
  Twf_cur = new double[N];
  Twf_new = new double[N];
  tmp     = new double[N];
  double *swaptmp;

  // Twf_old=f
  for(int n=0;n<N;n++)
    Twf_old[n]=f[n];
  // Twf_cur = (L*f-a2*f)/a1;
  L->sparse_matrix_multiply(f,Twf_cur);
  for(int n=0;n<N;n++)
    Twf_cur[n]=(Twf_cur[n]-a2*f[n])/a1;

  for(int ks=0;ks< c->s;ks++){
    //  r{j}=.5*c{j}(1)*Twf_old + c{j}(2)*Twf_cur;
    for(int n=0;n<N;n++)
      res[ks][n]=.5*c->c[ks][0]*Twf_old[n] + c->c[ks][1]*Twf_cur[n];
  }
  for(int km=1;km<c->m-1;km++){
    //  Twf_new = (2/a1)*(L*Twf_cur-a2*Twf_cur)-Twf_old;
    L->sparse_matrix_multiply(Twf_cur,tmp);
    for(int n=0;n<N;n++)
      Twf_new[n]=(2.0/a1)*(tmp[n]-a2*Twf_cur[n])-Twf_old[n];

    for(int ks=0;ks< c->s;ks++){
      //    r{j}=r{j}+c{j}(k+1)*Twf_new;
      for(int n=0;n<N;n++)
	res[ks][n]+=c->c[ks][km+1]*Twf_new[n] ;
    }
    // Twf_old=Twf_cur;
    // Twf_cur=Twf_new;
    swaptmp=Twf_old;
    Twf_old=Twf_cur;
    Twf_cur=Twf_new;
    Twf_new=swaptmp;
  }
  delete[] Twf_old;
  delete[] Twf_cur;
  delete[] Twf_new;
  delete[] tmp;
}
// y[ks] are coefficients at scale ks
// so y[ks] must be valid for 0<=ks< c->s, y[ks] points to N doubles
// res points to N doubles
void sgwt_cheby_op_adjoint(double *res, double *y[],sparse_compressed_row *L, cheby_coeff *c, approximation_range *arange){
  // 
  double *tmp=new double[L->N];
  // I want to create separate cheby_coeff instances for each scale,
  // call sgwt_cheby_op to apply these polynomials to their corresponding sgwt coefficient vectors, then accumulate
  cheby_coeff *c_single[c->s];
  for (int ks=0;ks<c->s;ks++)
    c_single[ks]=new cheby_coeff(c,ks);

  for (int n=0;n<L->N;n++)
    res[n]=0;

  for (int ks=0;ks<c->s;ks++){
    cheby_coeff c_one(c,ks);
    sgwt_cheby_op(&tmp,y[ks],L,c_single[ks],arange);
    for (int n=0;n<L->N;n++)
      res[n]+=tmp[n];
  }
 for (int ks=0;ks<c->s;ks++)
   delete c_single[ks];
 delete [] tmp;
}

bool cheby_op::is_valid_input(const mxArray *L_in,const mxArray *c_in,const mxArray *arange_in){
  return (sparse_compressed_row::is_valid_input(L_in) &&
	  cheby_coeff::is_valid_input(c_in) &&
	  approximation_range::is_valid_input(arange_in) ); 
}

cheby_op::cheby_op(const mxArray *L_in,const mxArray *c_in,const mxArray *arange_in){
  L=new sparse_compressed_row(L_in);
  c=new cheby_coeff(c_in);
  arange=new approximation_range(arange_in);
  N=L->N;
  s=c->s;
  // allocate
  Twf_old=new double[N];
  Twf_cur=new double[N];
  Twf_new=new double[N];
  Twf_tmp=new double[N];

  tmp_res=new double[N];
  c_single.resize(c->s);
  for(int ks=0;ks<c->s;ks++)
    c_single[ks]=new cheby_coeff(c,ks);
}

cheby_op::~cheby_op(){
  delete L;
  delete c;
  delete arange;

  delete [] Twf_old;
  delete [] Twf_cur;
  delete [] Twf_new;
  delete [] Twf_tmp;
  delete [] tmp_res;  
  for(unsigned int ks=0;ks < c_single.size();ks++)
    delete c_single[ks];
}

void cheby_op::cheby_op_forward(double *res[],double *f){
  double a1=((*arange)[1]-(*arange)[0])/2;
  double a2=((*arange)[1]+(*arange)[0])/2;
  double *swaptmp;
  /* 
     This code is hard to read due to low-level c-style array syntax.
     It is computing the chebyshev polynomials of L, applied to f,
     and storing these in res

     res[ks] is N-length vector of doubles for result at scale ks
     c->c[ks][km+1] is 1+km degree Chebyshev coefficient for the ks'th polynomial.
  */
  // Twf_old=f
  for(int n=0;n<N;n++)
    Twf_old[n]=f[n];
  // Twf_cur = (L*f-a2*f)/a1;
  L->sparse_matrix_multiply(f,Twf_cur);
  for(int n=0;n<N;n++)
    Twf_cur[n]=(Twf_cur[n]-a2*f[n])/a1;

  for(int ks=0;ks< c->s;ks++){
    //  r{j}=.5*c{j}(1)*Twf_old + c{j}(2)*Twf_cur;
    for(int n=0;n<N;n++)
      res[ks][n]=.5*c->c[ks][0]*Twf_old[n] + c->c[ks][1]*Twf_cur[n];
  }
  for(int km=1;km<c->m-1;km++){
    //  Twf_new = (2/a1)*(L*Twf_cur-a2*Twf_cur)-Twf_old;
    L->sparse_matrix_multiply(Twf_cur,Twf_tmp);
    for(int n=0;n<N;n++)
      Twf_new[n]=(2.0/a1)*(Twf_tmp[n]-a2*Twf_cur[n])-Twf_old[n];

    for(int ks=0;ks< c->s;ks++){
      //    r{j}=r{j}+c{j}(k+1)*Twf_new;
      for(int n=0;n<N;n++)
	res[ks][n]+=c->c[ks][km+1]*Twf_new[n] ;
    }
    // Twf_old=Twf_cur;
    // Twf_cur=Twf_new;
    swaptmp=Twf_old;
    Twf_old=Twf_cur;
    Twf_cur=Twf_new;
    Twf_new=swaptmp;
  }
}

void cheby_op::cheby_op_adjoint(double *res,double *y[]){
  for (int n=0 ; n <N;n++ )
    res[n]=0;
  // Apply each single chebyshev polynomial of L to corresponding 
  // coefficients, and accumulate result
  for (int ks=0;ks<c->s;ks++){
    sgwt_cheby_op(&tmp_res,y[ks],L,c_single[ks],arange);
    for (int n=0;n<L->N;n++)
      res[n]+=tmp_res[n];
  }
}
