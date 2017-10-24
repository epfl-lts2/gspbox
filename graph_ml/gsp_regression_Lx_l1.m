function sol = gsp_regression_Lx_l1(G, M, y, tau, param )
%GSP_REGRESSION_LX_L1 Regression using a piecewise linear assumption
%   Usage: sol = gsp_regression_tv(G, M, y, tau );
%          sol = gsp_regression_tv(G, M, y, tau, param );
%
%   Input parameters:
%       G   : Graph
%       M   : Mask (to determine which label is known)
%       y   : label (total size of the problem)
%       tau : regularization parameter (weight for tv)
%       param : optional structure of parameters
%
%   Output parameters:
%       sol : Solution of the problem
%
%   This function solve the following problem
%
%   .. argmin_x  || M x - y ||_2^2 + tau || Lx ||_1
%
%   If tau is set to zero, then the following problem is solved
%
%   ..  argmin_x   || Lx ||_1   s. t.  M x - y = 0
%
%   This function uses the UNLocBoX.
%
%   See also: gsp_regression_tv gsp_regression_tik gsp_classification_tv

% Author: Vassilis Kalofolias
% Date  : february 2016


%% Optional parameters

if nargin < 5
    param = struct;
end

if nargin < 4
    tau = 0;
end


if ~isfield(param,'verbose'), param.verbose = 1; end

%% prepare the graph

G = gsp_estimate_lmax(G);
G = gsp_adj2vec(G);





paramsolver = param;

    
if tau > 0
    M_op =@(x) bsxfun(@times, M, x);

    fg.grad = @(x) 2 * M_op(M_op(x) - y);
    fg.eval = @(x) norm(M_op(x) - y)^2;
    fg.beta = 2;
    
    % setting the function ||Lx||_1
    paramtv.verbose = param.verbose-1;
    ftv.prox = @(x,T) prox_l1(x,tau*T,paramtv);
    ftv.eval = @(x) tau *sum(vec(G.L*x)); 
    ftv.L = @(x ) G.L * x;
    ftv.Lt = @(x ) G.L * x;
    ftv.norm_L = G.lmax^2;
    %% solve the problem
    % setting the timestep

    
    sol = solvep(y, {ftv, fg}, paramsolver);


else
% %     param_b2.verbose = param.verbose -1;
% %     param_b2.y = y;
% %     param_b2.A = @(x) M.*x;
% %     param_b2.At = @(x) M.*x;
% %     param_b2.tight = param.tight;
% %     param_b2.epsilon = 0;
% %     fproj.prox = @(x,T) proj_b2(x,T,param_b2);
% %     fproj.eval = @(x) eps;
% 
%     fproj.prox = @(x,T) x - M.*x + M.*y;
%     fproj.eval = @(x) eps;
%     
%     % setting the function ftv
% 
% %     paramtv.verbose = param.verbose-1;
% %     ftv.prox = @(x,T) gsp_prox_tv(x,T,G,paramtv);
% %     ftv.eval = @(x) sum(gsp_norm_tv(G,x));  
%     
%     paramtv.verbose = param.verbose-1;
%     ftv.prox = @(x,T) prox_l1(x,T,paramtv);
%     ftv.eval = @(x) sum(gsp_norm_tv(G,x)); 
%     ftv.L = @(x ) G.Diff * x;
%     ftv.Lt = @(x ) G.Diff' * x;
%     ftv.norm_L = G.lmax;
%    
%     %% solve the problem
%     sol = solvep(y,{ftv,fproj},paramsolver);


    A = G.L(:,~logical(M));
    b = - G.L(:,logical(M)) * y(logical(M),:);
    
    paramtv.verbose = param.verbose-1;
    paramtv.y = b;
    ftv.prox = @(x,T) prox_l1(x, T, paramtv);
    ftv.eval = @(x) sum(sum(abs(A*x-b))); 
    ftv.L = @(x) A * x;
    ftv.Lt = @(x) A' * x;
    ftv.norm_L = G.lmax^2;
    
    
    %% solve the problem
    solt = solvep(y(~logical(M),:), {ftv}, paramsolver);
    
    sol = y;
    sol(~logical(M), :) = solt;

end






% sol = sol(logical(1-M));
% sol = reshape(sol,[],size(M,2));

end
