function r = gsp_jtv_cheby_op(G, c, signal,param)
%GSP_JTV_CHEBY_OP : Chebyshev polynomial of graph Laplacian applied to time-vertex signals
%   Usage: r = gsp_jtv_cheby_op(G, c, signal)
%              gsp_jtv_cheby_op(G, c, signal,param)
%
%   Input parameters:
%       G       : Time-Vertex Graph structure
%       c       : Chebyshef coefficients
%       signal  : Time-Vertex signal
%       param   : Structure of optional parameters
%   Output parameters
%       r       : Coefficients matrix
% 
%
%   This function is inspired by the function gsp_cheby_op
%
%   See also: gsp_jtv_cheby_coeff, gsp_jtv_filter_analysis
%

% Author: Francesco Grassi
% Date : July 2016
% Testing: test_jtv_filter

if nargin < 4
    param = struct;
end

if ~isfield(param,'verbose'), param.verbose = 1; end;


Nscales=size(c,3);

M = size(c,1);


if ~isfield(G,'lmax');
    G = gsp_estimate_lmax(G);
    if param.verbose
    warning(['GSP_JTV_CHEBY_OP: The variable lmax is not ',...
        'available. The function will compute it for you. ',...
        'However, if you apply many time this function, you ',...
        'should precompute it using the function: ',...
        'gsp_estimate_lmax']);
    end
end


arange = [0, G.lmax];

a1 = (arange(2) - arange(1))/2;
a2 = (arange(2) + arange(1))/2;

T = G.jtv.T;

if G.jtv.extension
    tau = 2*T-1;
else
    tau = T;
end

Twf_old=signal;                    
Twf_cur=(G.L*Twf_old-a2*Twf_old)/a1; 

r = zeros(G.N,tau,Nscales);

for ii=1:Nscales
    r(:,:,ii) = 0.5 * Twf_old.*repmat(c(1,:,ii),G.N,1) + Twf_cur.*repmat(c(2,:,ii),G.N,1);
end

for k=2:M
    Twf_new = (2/a1) * (G.L*Twf_cur-a2*Twf_cur) - Twf_old;
    for ii=1:Nscales
        if 1+k <= M
            r(:,:,ii) = r(:,:,ii) + Twf_new.*repmat(c(k+1,:,ii),G.N,1);
        end
    end
    Twf_old=Twf_cur;
    Twf_cur=Twf_new;
end


end
