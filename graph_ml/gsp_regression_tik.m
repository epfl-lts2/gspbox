function sol = gsp_regression_tik(G ,M, y , tau, param )
%GSP_REGRESSION_TIK Regression using graph and Tikhonov
%   Usage: sol = gsp_regression_tik(G ,M, y , tau );
%          sol = gsp_regression_tik(G ,M, y , tau, param );
%
%   Input parameters:
%       G   : Graph
%       M   : Mask (to determine with label is known)
%       y   : label (total size of the problem)
%       tau : regularization parameter (weight for tv)
%       param : optional structure of parameters
%
%   Output parameters:
%       sol : Solution of the problem
%
%   This function solve the following problem
%
%   .. argmin_x  || M x - y ||_2^2 + tau || nabla_G x ||_2^2
%
%   If tau is set to zero, then the following problem is solved
%
%   ..  argmin_x   || nabla_G x ||_2^2    s. t.  M x - y = 0
%
%   For the las problem, this function can compute an exact solution if
%   *param.exact* is activated. It will be efficient if the number of
%   unlabelled points is low.
%
%   Additional parameters
%   ---------------------
%
%   * *param.verbose* : Verbosity of the iterative algorithm
%   * *param.maxit* : maximum number of iteration for PCG
%
%   This function uses the UNLocBoX. 
%
%   See also: gsp_classification_tik gsp_regression_tv
%

% Author: Nathanael Perraudin
% Date  : 24 July 2015
% Testing: test_graph_ml



%% Optional parameters

if nargin<5
    param = struct;
end

if nargin<4
    tau = 0;
end



if ~isfield(param,'verbose'), param.verbose = 1; end
if ~isfield(param,'tol'), param.tol = 1e-6; end
if ~isfield(param,'maxit'), param.maxit = 200; end
% if ~isfield(param,'direct'), param.direct = (tau==0); end
% if ~isfield(param,'exact'), param.exact = (numel(M)-nnz(M))<1000; end

if tau==0
    if param.verbose
        fprintf('Using direct solution \n')
    end
    if (numel(M) == size(M,1)) || (numel(M) == size(M,2))
        indl = find(M);
        indu = find(1-M);        
    else   
        error('I cannot handle this case yet');
    end
    Luu = (G.L(indu,indu));
    Wul = - G.L(indu,indl);
    tmp = (Wul * y(indl,:));
%     if ~param.exact
%         if ~isfield(param,'order'), param.order = 30; end
%         Gtemp.L = Luu;
%         Gtemp.N = size(Luu,1);
% %        Gtemp = gsp_estimate_lmax(Gtemp);
% %         cheb_coeffs = gsp_cheby_coeff(Gtemp, @(x) pinv_n(x,1e-8),...
% %         param.order, param.order +1);
% %         solt = gsp_cheby_op(Gtemp, cheb_coeffs, tmp);
%         paramt.method = 'lanczos';
%         paramt.order = param.order;
%         solt = gsp_filter_analysis(Gtemp, @(x) pinv_n(x,1e-8),tmp,paramt);
% 
%     else
%        solt = pinv(full(Luu)) * tmp;
         solt = Luu \ tmp;
%     end
    sol = y;
    sol(indu,:) = solt;
    return
else


    %% prepare the graph
    if ~isfield(G,'lmax')
        G = gsp_estimate_lmax(G);
    end

    %% set the 
    % setting the function f2 (see unlocbox for help)

    Mop =@(x) bsxfun(@times,M,x);

%     fg.grad = @(x) 2*Mop(Mop(x)-y);
%     fg.eval = @(x) norm(Mop(x)-y)^2;
%     fg.beta = 2;
% %     paramtik.verbose = param.verbose -1;
% %     ftik.prox = @(x,T) gsp_prox_tik(x,tau * T,G,paramtik);
%     ftik.eval = @(x) tau* sum(gsp_norm_tik(G,x));
%     ftik.grad = @(x) 2 * G.L * x;
%     ftik.beta = 2 * G.lmax;
    
% 
% 
% else
% %     param_b2.verbose = param.verbose -1;
% %     param_b2.y = y;
% %     param_b2.A = @(x) M.*x;
% %     param_b2.At = @(x) M.*x;
% %     param_b2.tight = param.tight;
% %     param_b2.epsilon = 0;
% %     fproj.prox = @(x,T) proj_b2(x,T,param_b2);
% %     fproj.eval = @(x) eps;
% 
%     fproj.prox = @(x,T) x - Mop(x) + Mop(y);
%     fproj.eval = @(x) eps;
%     ftik.eval = @(x) sum(gsp_norm_tik(G,x));   
%     Ltmp = G.L + G.L';
%     ftik.grad = @(x) Ltmp*x;
%     ftik.beta = 2*G.lmax;
% end

% %% solve the problem
% 
% % setting different parameter for the simulation
% paramsolver = param;
% 
% if tau > 0
%     sol = forward_backward(y,ftik,fg,paramsolver);
% else
%     sol = forward_backward(y,fproj,ftik,paramsolver);
% end

% sol = sol(logical(1-M));
% 
% sol = reshape(sol,[],size(M,2));

    A = @(z)  vec ( Mop(reshape(z,[],size(y,2))) ...
        +  tau * G.L * reshape(z,[],size(y,2)));
    b =  vec(Mop(y));        
    sol = pcg(A,b,param.tol,param.maxit,[],[],b);
    sol = reshape(sol,[],size(y,2));
end

end


% function r =  pinv_n(x,t)
% 
% r = double(abs(x)>t) .* 1./x;
% 
% end


%   * *param.direct* : Direct computation of the exact solution (only for
%     $\tau = 0$). (Default tau==0)
%   * *param.exact* : Exact computation of the exact solution (only for
%     $\tau = 0$ and param.direct = 0). (Default: (numel(M)-nnz(M))<1000 )