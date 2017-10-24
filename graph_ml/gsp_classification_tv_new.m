function sol = gsp_classification_tv_new(G ,M, y , tau, param )
%GSP_CLASSIFICATION_TV Classification using graph and Tikonow
%   Usage: sol = gsp_classification_tv(G ,M, y );
%          sol = gsp_classification_tv(G ,M, y , tau );
%          sol = gsp_classification_tv(G ,M, y , tau, param );
%
%   Input parameters:
%       G   : Graph
%       M   : Mask (to determine with label is known)
%       y   : label (total size of the problem)
%       tau : regularization parameter (weight for tv) (default 0)
%       param : optional structure of parameters
%
%
%   Output parameters:
%       sol : Solution of the problem
%
%   This function solve the following problem
%
%   .. argmin_x  || M x - B ||_2^2 + tau || x ||_{G TV}
%
%   where B is a matrix create using the function |gsp_classification_matrix|  
%
%   If tau is set to zero, then the following problem is solved:
%
%   ..  argmin_x   || x ||_{G TV}    s. t.  M x - B = 0
%   
%   Warning the class needs to be integers! (Consecutive for optimality)
%
%   This function uses the UNLocBoX.
%
%   See also: gsp_regression_tv gsp_classification_tik

% Author: Nathanael Perraudin
% Date  : 24 July 2015



%% Optional parameters

if nargin<5
    param = struct;
end

if nargin<4
    tau = 0;
end


if ~isfield(param,'verbose'), param.verbose = 1; end

%%
B = gsp_classification_matrix(y);


%% prepare the graph

G = gsp_estimate_lmax(G);
G = gsp_adj2vec(G);





paramsolver = param;

    
if tau > 0
%     Mop =@(x) bsxfun(@times,M,x);
% 
%     fg.grad = @(x) 2*Mop(Mop(x)-y);
%     fg.eval = @(x) norm(Mop(x)-y)^2;
%     fg.beta = 2;
%     
%     % setting the function ftv
% %     paramtv.verbose = param.verbose-1;
% %     ftv.prox = @(x,T) gsp_prox_tv(x,tau*T,G,paramtv);
% %     ftv.eval = @(x) tau *sum(gsp_norm_tv(G,x));   
% 
%     paramtv.verbose = param.verbose-1;
%     ftv.prox = @(x,T) prox_l1(x,tau*T,paramtv);
%     ftv.eval = @(x) tau *sum(gsp_norm_tv(G,x)); 
%     ftv.L = @(x ) G.Diff * x;
%     ftv.Lt = @(x ) G.Diff' * x;
%     ftv.norm_L = G.lmax;
%     %% solve the problem
%     % setting the timestep
% 
%     
%     sol = solvep(y,{ftv,fg},paramsolver);

error('not done')

else


    A = G.Diff(:,~logical(M));
    b = - G.Diff(:,logical(M)) * B(logical(M),:);
    
    paramtv.verbose = param.verbose-1;
    paramtv.y = b;
    ftv.prox = @(x,T) prox_l1(x,T,paramtv);
    ftv.eval = @(x) sum(sum(abs(A*x-b))); 
    ftv.L = @(x ) A * x;
    ftv.Lt = @(x ) A' * x;
    ftv.norm_L = G.lmax;    
    
    paramsimplex.dim = 2;
    paramsimplex.verbose = param.verbose - 1;
    fsimplex.prox = @(x,T) proj_simplex(x,T, paramsimplex);
    fsimplex.eval = @(x) eps;
    
    %% solve the problem
    solt = solvep(B(~logical(M),:),{ftv, fsimplex},paramsolver);
    
    soltv = B;
    soltv(~logical(M),:) = solt;

end


minf = min(y);
sol = gsp_matrix2label(soltv,minf);



end