function [G]=gsp_spectrum_cdf_approx(G,param)
%GSP_SPECTRUM_CDF_APPROX  Compute an approximation of the cumulative density function of the graph Laplacian eigenvalues
%   Usage:  G=gsp_spectrum_cdf_approx(G);
%           G=gsp_spectrum_cdf_approx(G,param);
%
%   Input parameters:
%         G                     : Graph structure.
%   Output parameters:
%         G                     : Graph structure, including the addition of G.spectral_warp_fn
%         param                 : Structure of additional parameter
%   Additional parameters:
%         param.num_pts         : Number of interpolation points
%         param.use_speedup     : Only perform step 1 (ldlsymbol) of the ldl algorithm once and then repeat the numeric step 2
%         param.use_permutation : Use a sparsity-preserving permutation 
%         param.use_ldl_package : If the mex files for Timothy Davis' ldlsparse package are present. Otherwise use MATLAB's 'ldl' function
%         param.ldl_thresh      : Threshold parameter for MATLAB 'ldl' function
%
%
%   'gsp_spectrum_cdf_approx(G)' computes an approximation of the
%   cumulative density function of the graph Laplacian eigenvalues using 
%   spectrum slicing techniques. The algorithm is as follows:
%
%   Step 1: Take $Q+1$ (param.num_pts) evenly spaced points on the interval
%   $[0,\lambda_{\max}]$.
%
%   Step 2: For every $q \in \{1,2,\ldots,Q\}$, compute a triangular
%   factorization:
%
%   .. \{\cal L}-\frac{q \lambda_{\max}}{Q}I = L_q \Delta_q L_q^*
%
%   .. math:: {\cal L}-\frac{q \lambda_{\max}}{Q}I = L_q \Delta_q L_q^* 
%
%   where $\Delta_q$ is a diagonal matrix and $L_q$ is a lower triangular
%   matrix. By a corollary of Sylverster's law of inertia, the number of
%   negative eigenvalues of the diagonal matrix $\Delta_q$ is equal to the
%   number of eigenvalues of ${\cal L}$ less than $q \lambda_{\max} / Q$.
%
%   Step 3: Use monotonic cubic polynomial interpolation to interpolate the
%   points
%
%   .. { ( q\lambda_{max}/Q , mu_q/(N-1) ) }_q=0,1,\ldots,Q
%
%   .. math:: \left\{\left(\frac{q \lambda_{upper}}{Q},\frac{\mu_q}{N-1} \right)\right\}_{q=0,1,\ldots,Q}
%
%   where $\mu_q$ is the number of diagonal elements of $\Delta_q$ less than zero. 
%
%   To perform the triangular factorizations, we use the LDL sparse
%   Cholesky package, written by Timothy Davis. The option
%   param.use_permutation enables a sparsity preserving permutation with
%   the function 'symamd'. If the option param.use_speedup is set to 1, the
%   symbolic analysis step of the LDL algorithm is only performed once for
%   all shifts.
%
%   See also:  
%
%   Demos:  
%
%   Requires: ldlsymbol_extra, ldlnumeric, ldlsparse
% 
%   References: parlett1998symmetric davis2005ldl davis2011ldl
%

%TODO this function need to be cleaned

%   AUTHOR : David I Shuman.
%   TESTING: 
%   REFERENCE:
  
if isfield(G, 'spectral_warp_fn')
    warning('Overwriting spectral warping function');
end

% if isfield(G,'e');
%     [x,y] = gsp_point2dcdf(G.e);
%     G.spectrum_cdf_approx = @(s) gsp_mono_cubic_warp_fn(x,y,s);
%     return
% end

if nargin<2
    param = struct;
end

if ~isfield(param, 'num_pts')
    num_pts=25;
else
    num_pts = param.num_pts;
end

if ~isfield(param, 'use_speedup')
    if ( exist('ldlsymbol_extra','file')==3 && exist('ldlnumeric','file')==3 )
        use_speedup=1;
    else
        use_speedup=0;
    end
else
    use_speedup = param.use_speedup;
end

if ~isfield(param, 'use_ldl_package')
    if exist('ldlsparse','file')==3
        use_ldl_package=1;
    else
        use_ldl_package=0;
    end
else
    use_ldl_package = param.use_ldl_package;
end

if ~isfield(param, 'use_permutation')
    use_permutation=1;
else
    use_permutation = param.use_permutation;
end

if ~isfield(G, 'lmax')
    warning(['GSP_SPECTRUM_CDF_APPROX: The variable lmax is not ',...
            'available. The function will compute it for you. ',...
            'However, if you apply many time this function, you ',...
            'should precompute it using the function: ',...
            'gsp_estimate_lmax']);
    G = gsp_estimate_lmax(G);
end

if ~isfield(param, 'ldl_thresh')
    ldl_thresh=.001;
else
    ldl_thresh = param.ldl_thresh;
end



counts=zeros(num_pts,1);
counts(num_pts)=G.N-1;

interp_x=(0:num_pts-1)*G.lmax/(num_pts-1);
interp_x=interp_x';

identity=speye(G.N);

if use_speedup
    if use_permutation
        P=symamd(G.L);
        [Parent, Lp, PO, PIn, flopcount] = ldlsymbol_extra(G.L,P);
        for i=1:num_pts-2 
            mat=G.L-interp_x(i+1)*identity;
            [~, HD]=ldlnumeric(mat,Lp,Parent,PO,PIn);
            counts(i+1)=sum(diag(HD)<0);
        end
    else
        [Parent, Lp, flopcount] = ldlsymbol_extra(G.L);
        for i=1:num_pts-2 
            mat=G.L-interp_x(i+1)*identity;
            [~, HD]=ldlnumeric(mat,Lp,Parent);            
            counts(i+1)=sum(diag(HD)<0);
        end
    end
else
    if use_permutation
        P=symamd(G.L);
        for i=2:num_pts-1 
            mat=G.L-interp_x(i)*identity;
            if use_ldl_package
                [~,HD]=ldlsparse(mat,P);
            else
                [~,HD,~]=ldl(mat(P,P),ldl_thresh);
            end
            counts(i)=sum(diag(HD)<0);
        end
    else
        for i=2:num_pts-1 
            mat=G.L-interp_x(i)*identity;
            if use_ldl_package
                [~,HD]=ldlsparse(mat);
            else
                [~,HD,~]=ldl(mat,ldl_thresh);
            end
            counts(i)=sum(diag(HD)<0);
        end
    end
end
interp_y=counts/(G.N-1);

G.spectrum_cdf_approx = @(s) gsp_mono_cubic_warp_fn(interp_x,interp_y,s);

end
