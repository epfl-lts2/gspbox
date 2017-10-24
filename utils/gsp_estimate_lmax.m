function G = gsp_estimate_lmax(G)
%GSP_ESTIMATE_LMAX estimates the maximum laplacian eigenvalue
%   Usage: G = gsp_estimate_lmax(G);
%
%   Inputs parameters:
%       G   : Graph structure (or cell array of graph structure)
%   Outputs parameters:
%       G   : Graph structure (or cell array of graph structure)
%
%   This function will compute an approximation of the maximum laplacian
%   eigenvalue and fill it in the field *G.lmax*
%
%   Runs Arnoldi algorithm with a large tolerance, then increases
%   calculated maximum eigenvalue by 1 percent. For much of the gspbox
%   machinery, we need to approximate the wavelet kernels on an
%   interval that contains the spectrum of L. The only cost of using
%   a larger interval is that the polynomial approximation over the
%   larger interval may be a slightly worse approxmation on the
%   actual spectrum. As this is a very mild effect, it is not likely
%   necessary to obtain very tight bonds on the spectrum of L
%
%   This function is inspired by the sgwt_toolbox

% Author: David K. Hammond, Nathanael Perraudin
% Date:   15 March 2014


if numel(G)>1
    Ng = numel(G);
    for ii = 1:Ng
       G{ii} = gsp_estimate_lmax(G{ii});
    end     
    return;
end


try
    opts=struct('tol',5e-3,'p',min(G.N,10),'disp',0);
    lmax=eigs(G.L,1,'lm',opts);

    G.lmax=abs(lmax)*1.01; % just increase by 1 percent to be robust to error
    
catch
    warning('GSP_ESTIMATE_LMAX: Cannot use the default method')
    G.lmax = 2*max(G.d);
end

if isfield(G,'Gm')
    G = gsp_estimate_oose_lmax(G);
end

end