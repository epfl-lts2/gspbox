function r = gsp_cheby_op(G, c, signal,param)
%GSP_CHEBY_OP : Chebyshev polynomial of graph Laplacian applied to vector
%   Usage: r = gsp_cheby_op(G, c, signal)
%
%   Input parameters:
%       G       : Graph structure
%       c       : Chebyshef coefficients
%       signal  : Signal to filter
%   Output parameters
%       r       : Result of the filtering
% 
%   Compute (possibly multiple) polynomials of graph laplacian (in
%   Chebyshev basis) applied to input.
%
%   Coefficients for multiple polynomials may be passed as a matrix.
%   This is equivalent to setting::
%
%       r(1) = gsp_cheby_op(G, c(:,1), signal);
%       r(2) = gsp_cheby_op(G, c(:,2), signal);
%       ...
% 
%   but is more efficient as the Chebyshev polynomials of G.L applied
%   to *signal* can be computed once and shared.
%
%   The output *r* is a matrix with each column corresponding to a filter.
%
%   *param* contain only one field param.verbose to controle the verbosity.
%
%   Example:::
%
%         Nf = 4;
%         G = gsp_sensor(100);
%         G = gsp_estimate_lmax(G);
%         g = gsp_design_meyer(G, Nf);  
%         c = gsp_cheby_coeff(G, g);
%         f = rand(G.N,1);
%         r = gsp_cheby_op(G, c, f);
%
%   This function is inspired by the sgwt_toolbox
%
%   See also: gsp_cheby_coeff gsp_filter_analysis
%

% Author: David K Hammond, Nathanael Perraudin
% Testing: test_filter
% Date: 19 March 2014

if nargin < 4
    param = struct;
end

if ~isfield(param,'verbose'), param.verbose = 1; end;


Nscales=size(c,2);

M = size(c,1);
% To handle different order of Cheby approximation


assert(all(M>=2));

maxM=max(M);



if ~isfield(G,'lmax');
    G = gsp_estimate_lmax(G);
    if param.verbose
    warning(['GSP_CHEBY_OP: The variable lmax is not ',...
        'available. The function will compute it for you. ',...
        'However, if you apply many time this function, you ',...
        'should precompute it using the function: ',...
        'gsp_estimate_lmax']);
    end
end

if isa(signal,'single')
    signal = double(signal); 
    bsingle = 1;
else
    bsingle = 0;
end

arange = [0, G.lmax];

a1 = (arange(2) - arange(1))/2;
a2 = (arange(2) + arange(1))/2;


%Twf_new = T_j(L) f
%Twf_cur T_{j-1}(L) f
%TWf_old T_{j-2}(L) f

Twf_old=signal;                     % j = 0;
Twf_cur=(G.L*signal-a2*signal)/a1;  % j = 1;

Nv = size(signal,2);
r = zeros(G.N*Nscales,Nv);

for ii=1:Nscales
    r((1:G.N)+G.N * (ii-1),:) = 0.5 * c(1,ii) * Twf_old + c(2,ii) * Twf_cur;
end

for k=2:maxM
    Twf_new = (2/a1) * (G.L*Twf_cur-a2*Twf_cur) - Twf_old;
    for ii=1:Nscales
        if 1+k <= M
            r((1:G.N)+G.N * (ii-1),:) =...
                r((1:G.N)+G.N * (ii-1),:) + c(k+1,ii)*Twf_new;
        end
    end
    Twf_old=Twf_cur;
    Twf_cur=Twf_new;
end


if bsingle
    r = single(r);
end

end
