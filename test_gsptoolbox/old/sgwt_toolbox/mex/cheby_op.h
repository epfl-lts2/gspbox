#ifndef CHEBY_OP_H
#define CHEBY_OP_H
#include "mex.h"
#include <vector>
using std::vector;
// forward declarations
class sparse_compressed_row;
class cheby_coeff;
class approximation_range;
/*
  cheby_op : class to implement forward and adjoint SGWT transformation

  Basic use (for now, in the context of mex programming) is to construct the class with
  matlab arrays for L (the graph Laplacian), c (cell array of chebyshev coefficients) and
  arange (the interval containing spectrum of L over which the chebyshev polynomial approximation 
  is valid).

  Once initialized, forward and adjoint transforms may be computed by calling
  cheby_op_forward and cheby_op_adjoint

*/  
class cheby_op{
 public:
  int N; // dimension of graph
  int s; // # of scales (including scaling function). Equivalently, the # of distinct chebyshev polynomials
  /* cheby_op_forward : compute forward transform.
     res: array of pointers s.t. res[ks] points to N elements, for each 0<=ks<s
     f : length N array
  */
  void cheby_op_forward(double *res[],double *f);

  /* cheby_op_adjoint : compute adjoint transform.
     y: array of pointers s.t. y[ks] points to N elements, for each 0<=ks<s
     res : length N array
  */  
  void cheby_op_adjoint(double *res,double *y[]);
  cheby_op(const mxArray *L_in,const mxArray *c_in,const mxArray *arange_in);

  /* is_valid_input : verifies that inputs are valid to construct the corresponding
     classes representing L,c and arange */
  static bool is_valid_input(const mxArray *L_in,const mxArray *c_in,const mxArray *arange_in);
  ~cheby_op();
 private:
  sparse_compressed_row *L;
  cheby_coeff *c;
  approximation_range *arange;
  //  used in cheby_op (these are all length N arrays)
  double *Twf_old ;
  double *Twf_cur ;
  double *Twf_new ;
  double *Twf_tmp ; 
  // used in cheby_op_adjoint
  double *tmp_res ; // (length N array)
  vector<cheby_coeff *> c_single; // will hold cheby_coeff instances for single polynomials
};


/*
  sparse_compressed_row : class for representing sparse matrix in compressed-row format

  This essentially uses same representation as MATLAB's sparse matrix, but with the following difference:
  MATLAB uses compressed-column format. The constructor of this class calls MATLAB to take matrix transpose,
  so that the resulting pr,ir,jc data can be interpreted as compressed row format for original matrix.
*/

class sparse_compressed_row{
 public:
  int M;
  int N;
  sparse_compressed_row(const mxArray *array_in);
  ~sparse_compressed_row();
  // sparse_matrix_multiply : multiply vector x on left by matrix, storing result in res
  // prx must be of size N, res must be of size M. 
  void sparse_matrix_multiply(const double *x,double *res);
  // is_valid_input : checks that array_in is MATLAB sparse double matrix type
  static bool is_valid_input(const mxArray *array_in);
 private:
  int nz;
  double *pr;
  mwIndex *ir;
  mwIndex *jc;
};

/* 
   cheby_coeff : class for representing sets of polynomials in Chebyshev basis
   
   All polynomials are of the same degree, m
*/

class cheby_coeff{
 public:
  int m; // polynomial degree
  int s; // # of scales (# of distinct chebychev polynomials)
  vector<vector<double> > c; // the coefficients themselves
  // c[ks][k] is degree k coefficient for ks'th polynomial
  cheby_coeff(const mxArray *array_in);
  // this constructor creates new cheby_coeff containing a single specified
  // polynomial (this is useful for computing adjoint)
  cheby_coeff(cheby_coeff *c_in,int ks);
  // is_valid_input : checks that array_in is cell array, and that all cells are
  // double arrays with same # of elements
  static bool is_valid_input(const mxArray *array_in);
};
/*
  approximation_range : class for representing approximation range for chebyshev polynomials
*/
class approximation_range{
 public:
  approximation_range(const mxArray *array_in);
  // no boundary checking done
  double operator[](int k);
  static bool is_valid_input(const mxArray *array_in);
 private:
  double arange_[2];
};



void sgwt_cheby_op(double *res[],double *f,sparse_compressed_row *L, cheby_coeff *c, approximation_range *arange);
void sgwt_cheby_op_adjoint(double *res,double *y[],sparse_compressed_row *L, cheby_coeff *c, approximation_range *arange);

#endif
