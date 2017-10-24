function [s] = gsp_filter_inverse(G, filter, c, param)
%GSP_FILTER_INVERSE Inverse operator of a gsp filterbank
%   Usage:  s = gsp_filter_inverse(G, filter, c);
%           s = gsp_filter_inverse(G, filter, c, param);
%
%   Input parameters:
%         G         : Graph structure.
%         filter    : Set of spectral graph filters.
%         c         : Transform coefficients
%         param     : Optional parameter
%   Output parameters:
%         signal    : sythesis signal
%
%   'gsp_filter_inverse(G,filters,c)' computes the inverse
%   operator for coefficients $c$, where the atoms of the transform 
%   dictionary are generalized translations of each graph spectral filter
%   to each vertex on the graph.
%
%   .. f = (D'*D)^(-1) * D * c 
%
%   .. math:: f =  (D^*D)^{-1} D c
%
%   where the columns of $D$ are $g_{i,m}=T_i g_m$, and $T_i$ is a
%   generalized translation operator applied to each filter 
%   $\hat{g}_m(\cdot)$.  
%
%   Each column of *c* is the response of the signal to one filter.
%
%   Example:::
%
%         Nf = 4;
%         G = gsp_sensor(30);
%         G = gsp_estimate_lmax(G);
%         G = gsp_estimate_lmax(G);
%         g = gsp_design_mexican_hat(G, Nf);  
%         f = rand(G.N,1);
%         f = f/norm(f);
%         ff = gsp_filter_analysis(G,g,f);
%         f2 = gsp_filter_inverse(G,g,ff);
%         norm(f-f2)     
%
%   For additional parameters, please see |gsp_filter_synthesis|
%
%   See also: gsp_filter_analysis gsp_filter_synthesis
% 
%   References: hammond2011wavelets
%

% Author: Nathanael Perraudin
% Testing: test_filter
% Date: 19 March 2014

if nargin<4
    param = struct;
end

dual_filter = gsp_design_can_dual(filter);

s = gsp_filter_synthesis(G,dual_filter,c,param);


end



% 
% 
% %%% TODO: add lanczos and exact computation
% warning('This function solves the problem in a very stupid way')
% 
% % Read input parameters
% if nargin < 4
%     param = struct;
% end
% 
% Nf = length(filter);
% 
% if ~isfield(param,'order'); param.order = 30; end
% if ~isfield(param,'verbose'); param.verbose = 1; end
% if ~isfield(param,'tol'); param.tol = 1e-6; end
% if ~isfield(param,'maxit'); param.maxit = 200; end
% 
% if isfield(param,'method')
%     if strcmp(param.method,'lanczos');
%         error('Not implemented yet!');
%     end
% end
% 
%  
% if ~isfield(G,'lmax');
%     G = gsp_estimate_lmax(G);
%     if param.verbose
%         warning(['GSP_FILTER_ANALYSIS: The variable lmax is not ',...
%             'available. The function will compute it for you. ',...
%             'However, if you apply many time this function, you ',...
%             'should precompute it using the function: ',...
%             'gsp_estimate_lmax']);
%     end
% end
% 
% 
% cheb_coeffs = gsp_cheby_coeff(G, filter,...
%         param.order, param.order +1);    
% 
% 
% % Compute the adjoint
% adj = gsp_filter_synthesis(G, filter, c, param);
% 
% % W^* W
% % compute P(x) = p(x)^2
% M=size(cheb_coeffs,1);
% 
% d=zeros(1+2*(M-1),1);
% for ii=1:Nf
%     d=d+sgwt_cheby_square(cheb_coeffs(:,ii))';
% end
% wstarw = @(x) gsp_cheby_op(G,d,x);
% 
% % conjugate gradients
% s=zeros(G.N,size(adj,2));
% for ii=1:size(adj,2)
%     [s(:,ii),~,~,~]=pcg(wstarw,adj(:,ii),param.tol,param.maxit);
% end

