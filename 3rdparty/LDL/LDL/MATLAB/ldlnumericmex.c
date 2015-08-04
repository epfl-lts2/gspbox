/* ========================================================================== */
/* === ldlmex.c:  LDL mexFunction =========================================== */
/* ========================================================================== */

/* Code modified by David I Shuman to take the symbolic pattern as an input, 
 * in order to allow the user to just do the symbolic package once and then repeat the numeric part */

/* MATLAB interface for numerical LDL' factorization using the LDL sparse matrix
 * package.
 *
 * MATLAB calling syntax is:
 *
 *       [L, D] = ldlnumeric (A, Lp, Parent)
 *       [L, D] = ldlnumeric (A, Lp, Parent, P, PInv)
 *
 * The factorization is L*D*L' = A or L*D*L' = A(P,P).   A must be sparse,
 * square, and real.  L is lower triangular with unit diagonal, but the diagonal
 * is not returned.  D is diagonal sparse matrix.  Let n = size (A,1).   If P is
 * not present or empty, the factorization is:
 *
 *	(L + speye (n)) * D * (L + speye (n))' = A
 *
 * otherwise, the factorization is
 *
 *	(L + speye (n)) * D * (L + speye (n))' = A(P,P)
 *
 * P is a permutation of 1:n, an output of AMD, SYMAMD, or SYMRCM, for example.
 * Only the diagonal and upper triangular part of A or A(P,P) is accessed; the
 * lower triangular part is ignored.
 *
 * The elimination tree is returned in the Parent array.
 *
 * If no zero entry on the diagonal of D is encountered, then the flops argument
 * is the floating point count.
 *
 * If any entry on the diagonal of D is zero, then the LDL' factorization is
 * terminated at that point.  If there is no flops output argument, an error
 * message is printed and no outputs are returned.  Otherwise, flops is
 * negative, d = -flops, and D (d,d) is the first zero entry on the diagonal of
 * D.  A partial factorization is returned.  Let B = A if P is not present or
 * empty, or B = A(P,P) otherwise.  Then the factorization is
 *
 *	LDL = (L + speye (n)) * D * (L + speye (n))'
 *	LDL (1:d, 1:d) = B (1:d,1:d)
 *
 * That is, the LDL' factorization of B (1:d,1:d) is in the first d rows and
 * columns of L and D.  The rest of L and D are zero.
 *
 * Copyright (c) by Timothy A Davis, http://www.suitesparse.com.
 * All Rights Reserved.  See README for the License.
 */

#ifndef LDL_LONG
#define LDL_LONG
#endif

#include "ldl.h"
#include "mex.h"
#define Long SuiteSparse_long

/* ========================================================================== */
/* === LDL mexFunction ====================================================== */
/* ========================================================================== */

void mexFunction
(
    int	nargout,
    mxArray *pargout[ ],
    int	nargin,
    const mxArray *pargin[ ]
)
{
    Long i, n, *Pattern, *Flag, *Li, *Lp, *Ap, *Ai, *Lnz, *Parent,
	 lnz, do_num_only, *P, *Pinv, nn, k, j, permute, *Dp = NULL, *Di, d;
    double *Y, *D, *Lx, *Ax, *p, *lp ;

    /* ---------------------------------------------------------------------- */
    /* get inputs and allocate workspace */
    /* ---------------------------------------------------------------------- */

    do_num_only = ((nargin == 3) || (nargin == 5)) && (nargout <= 2) ;
    if (!do_num_only)
    {
	mexErrMsgTxt ("Usage:\n"
        "  [L, D] = ldlsparse (A, Lp, Parent) ;\n"
        "  [L, D] = ldlsparse (A, Lp, Parent, P, PInv)") ;
    }
    n = mxGetM (pargin [0]) ;
    if (!mxIsSparse (pargin [0]) || n != mxGetN (pargin [0])
	    || mxIsComplex (pargin [0]))
    {
    	mexErrMsgTxt ("ldl: A must be sparse, square, and real") ;
    }
    nn = (n == 0) ? 1 : n ;

    /* get sparse matrix A */
    Ap = (Long *) mxGetJc (pargin [0]) ;
    Ai = (Long *) mxGetIr (pargin [0]) ;
    Ax = mxGetPr (pargin [0]) ;
    
    /* get Lp */
    if (mxGetM (pargin [1]) * mxGetN (pargin [1]) != (n+1) ||
		mxIsSparse (pargin [1]))
	{
	    mexErrMsgTxt ("ldl: invalid input Lp\n") ;
	}
    Lp      = (Long *) mxMalloc ((n+1) * sizeof (Long)) ;
	lp = mxGetPr (pargin [1]) ;
	for (k = 0 ; k < (n+1) ; k++)
	{
	    Lp [k] = lp [k]; 
	}
    
    /* get Parent */
    if (mxGetM (pargin [2]) * mxGetN (pargin [2]) != n ||
		mxIsSparse (pargin [2]))
	{
	    mexErrMsgTxt ("ldl: invalid input Parent\n") ;
	}
    Parent  = (Long *) mxMalloc (nn * sizeof (Long)) ;
	p = mxGetPr (pargin [2]) ;
	for (k = 0 ; k < n ; k++)
	{
	    Parent [k] = p [k]; 
	}
        
    /* get fill-reducing ordering, if present */
    permute = ((nargin == 5) && !mxIsEmpty (pargin [3])) ;
    if (permute)
    {
	if (mxGetM (pargin [3]) * mxGetN (pargin [3]) != n ||
		mxIsSparse (pargin [3]))
	{
	    mexErrMsgTxt ("ldl: invalid input permutation\n") ;
	}
	P    = (Long *) mxMalloc (nn * sizeof (Long)) ;
	p = mxGetPr (pargin [3]) ;
	for (k = 0 ; k < n ; k++)
	{
	    P [k] = p [k]; /*p[k]-1 convert to 0-based */
	}
    if (mxGetM (pargin [4]) * mxGetN (pargin [4]) != n ||
        mxIsSparse (pargin [4]))
    {
        mexErrMsgTxt ("ldl: invalid input permutation inverse\n") ;
    }
    Pinv = (Long *) mxMalloc (nn * sizeof (Long)) ;
    p = mxGetPr (pargin [4]) ;
    for (k = 0 ; k < n ; k++)
    {
        Pinv [k] = p [k]; /*p[k]-1 convert to 0-based */
    }        
    }
    else
    {
	P    = (Long *) NULL ;
	Pinv = (Long *) NULL ;
    }

 
    /* get workspace */
    Lnz     = (Long *) mxMalloc (nn * sizeof (Long)) ;
    Flag    = (Long *) mxMalloc (nn * sizeof (Long)) ;
    Y       = (double *)  mxMalloc (nn * sizeof (double)) ;
    Pattern = (Long *) mxMalloc (nn * sizeof (Long)) ;
   
    /* make sure the input P is valid */
    if (permute && !ldl_l_valid_perm (n, P, Flag))
    {
	mexErrMsgTxt ("ldl: invalid input permutation\n") ;
    }

    /* note that we assume that the input matrix is valid */

    /* ---------------------------------------------------------------------- */
    /* symbolic factorization to get Lp, Parent, Lnz, and optionally Pinv */
    /* ---------------------------------------------------------------------- */
    /*ldl_l_symbolic (n, Ap, Ai, Lp, Parent, Lnz, Flag, P, Pinv) ;*/
    
    lnz = Lp [n] ;
    
    /* ---------------------------------------------------------------------- */
    /* create outputs */
    /* ---------------------------------------------------------------------- */

	/* create the output matrix L, using the Lp array from ldl_l_symbolic */
	pargout [0] = mxCreateSparse (n, n, lnz+1, mxREAL) ;
	mxFree (mxGetJc (pargout [0])) ;
	mxSetJc (pargout [0], (void *) Lp) ;	/* Lp is not mxFree'd */
	Li = (Long *) mxGetIr (pargout [0]) ;
	Lx = mxGetPr (pargout [0]) ;

	/* create sparse diagonal matrix D */
	if (nargout > 1)
	{
	    pargout [1] = mxCreateSparse (n, n, nn, mxREAL) ;
	    Dp = (Long *) mxGetJc (pargout [1]) ;
	    Di = (Long *) mxGetIr (pargout [1]) ;
	    for (j = 0 ; j < n ; j++)
	    {
		Dp [j] = j ;
		Di [j] = j ;
	    }
	    Dp [n] = n ;
	    D = mxGetPr (pargout [1])  ;
	}
	else
	{
	    D  = (double *) mxMalloc (nn * sizeof (double)) ;
	}
	
    /* ---------------------------------------------------------------------- */
    /* numeric factorization to get Li, Lx, and D */
    /* ---------------------------------------------------------------------- */
    
    d = ldl_l_numeric (n, Ap, Ai, Ax, Lp, Parent, Lnz, Li, Lx, D, Y, Flag,
	Pattern, P, Pinv) ;

    /* ---------------------------------------------------------------------- */
    /* singular case : truncate the factorization */
    /* ---------------------------------------------------------------------- */

    if (d != n)
    {
	    mexErrMsgTxt ("ldl: zero pivot encountered\n") ;
    }

    /* ---------------------------------------------------------------------- */
    /* free workspace */
    /* ---------------------------------------------------------------------- */

    if (nargout < 2)
    {
        mxFree (D) ;
    }
    if (permute)
    {
	mxFree (P) ;
	mxFree (Pinv) ;
    }
    mxFree (Parent) ;
    mxFree (Y) ;
    mxFree (Flag) ;
    mxFree (Pattern) ;
    mxFree (Lnz) ;
}
