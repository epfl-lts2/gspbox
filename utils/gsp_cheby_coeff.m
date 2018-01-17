function c = gsp_cheby_coeff(G, filter, m, N,param)
%GSP_CHEBY_COEFF : Compute Chebyshev coefficients for a filterbank
%   Usage: c = gsp_cheby_coeff(G, filter, m, N);
%          c = gsp_cheby_coeff(G, filter, m);
%          c = gsp_cheby_coeff(G, filter);
%
%   Input parameters:
%       G       : graph structure or range of application
%       filter  : filter or cell array of filters
%       m       : maximum order Chebyshev coefficient to compute (default 30)
%       N       : grid order used to compute quadrature (default is m+1)
%       param   : optional parameter
%   Output parameters
%       c   : matrix of Chebyshev coefficients
% 
%   This function compute the Chebyshef coefficients for all the filter
%   contained in the cell array filter. The coefficient are returned in a
%   matrix. Every collumn correspond to a filter. The coefficients are
%   ordered such that c(j+1) is j'th Chebyshev coefficient
%
%   Additional parameters
%   ---------------------
%  
%   * *param.use_chebfun*  : 1 to use the Chebfun package to compute 
%     Chebyshev coefficients
%   * *param.splitting_on* : 1 to call chebfun with splitting on
%   * *param.verbose* : Verbosity level (0 no log - 1 display warnings)
%     (default 1). 
%
%   Example:::
%
%         Nf = 4;
%         G = gsp_sensor(100);
%         G = gsp_estimate_lmax(G);
%         g = gsp_design_meyer(G, Nf);  
%         c = gsp_cheby_coeff(G, g);
%
%   This function is inspired by the sgwt_toolbox
%
%   See also: gsp_cheby_op gsp_filter_analysis
%

% Author: David K Hammond, Nathanael Perraudin, David Shuman
% Testing: test_filter
% Date: 19 March 2014

if nargin < 5
    param = struct;
end

if ~isfield(param,'verbose'), param.verbose = 1; end;

if nargin < 3
    m = 30;
end

if nargin < 4
   N = m+1; 
end

if iscell(filter)
   Nf = length(filter);
   c = zeros(m+1,Nf);
   for ii = 1: Nf
       c(:,ii) = gsp_cheby_coeff(G, filter{ii}, m, N,param);
   end
   return;
end

if isstruct(G)
    if ~isfield(G,'lmax');
        G = gsp_estimate_lmax(G);
        if param.verbose
        warning(['GSP_CHEBY_COEFF: The variable lmax is not ',...
            'available. The function will compute it for you. ',...
            'However, if you apply many time this function, you ',...
            'should precompute it using the function: ',...
            'gsp_estimate_lmax']);
        end
    end
  arange = [0, G.lmax];
else
  arange = G;
end

if ~isfield(param,'use_chebfun'), param.use_chebfun = 0; end;

if param.use_chebfun % Use Chebfun package, available at (http://www.chebfun.org/)
    if ~isfield(param,'splitting_on'), param.splitting_on = 0; end;
    if param.splitting_on
        h=chebfun(@(s) filter(s),arange,'splitting','on');
    else
        h=chebfun(@(s) filter(s),arange);
    end
    c=chebcoeffs(h,m+1); 
    c(1)=c(1)*2; 
else
    a1=(arange(2)-arange(1))/2;
    a2=(arange(2)+arange(1))/2;
    c = zeros(m+1,1);
    for ii=1:m+1
        c(ii) = sum( filter( a1* cos( (pi*((1:N)-0.5))/N) + a2) .* ...
                cos( pi*(ii-1)*((1:N)-.5)/N) ) *2/N;
    end
end

end
  
