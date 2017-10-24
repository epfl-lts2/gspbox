function [Lr, Sp, G, U, S, V] = gsp_frpcaog_1g(X, gamma, G, param)
%GSP_FRPCAOG_1G Fast robust PCA on graphs (1 graphs)
%   Usage: [Lr] = gsp_frpcaog_1g(X, gamma, G);
%          [Lr] = gsp_frpcaog_1g(X, gamma, G, param);
%          [Lr, Sp] = gsp_frpcaog_1g( ... );
%          [Lr, Sp, G] = gsp_frpcaog_1g( ... );
%          [Lr, Sp, G,  U, S, V] = gsp_frpcaog_1g( ... );
%
%   Input parameters: 
%       X       : Input data (matrix of double)
%       gamma   : Regularization parameter 1 (double)
%       G       : Graph (between the line of the data )
%       param   : Optional optimization parameters
%
%   Output Parameters:
%       Lr      : Low-rank part of the data
%       Sp      : Sparse part of the data
%       G       : Graph (between the line of the data )
%       U       : Part of the SVD of Lr
%       S       : Part of the SVD of Lr
%       V       : Part of the SVD of Lr
%
%   This function compute a low rank approximation of the data stored in
%   *Lr* by solving an optimization problem:
%
%   .. argmin_Lr || Lr - X ||_1  + gamma tr( X^T L X) 
%
%   .. math:: argmin_Lr || Lr - X ||_1  + \gamma tr( X^T L X) 
%
%   The sparse part is given by $ S = X - L_r $. 
%   
%   If $0$ is given for *G*, the corresponding graph will
%   be computed internally. The graph construction can be tuned using the
%   optional parameter: *param.paramnn*.
%
%   Other optional parameters are used in optimization. For details about
%   those, please see the help of the function |solvep| . 
%
%   If the number of output argument is greater than 2. The function, will
%   additionally compute a very economical SVD such that $ Lr = U S V^T$.
%
%   This function uses the UNLocBoX to be working.
%
%   References: shahid2015fast
%
%   See also: gsp_frpcaog_2g

% Author: Nathanael Perraudin
% Date  : 19th October 2015

%% Optional parameters

if nargin<4
    param = struct;
end

if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'paramnn'), param.paramnn = struct; end

if ~isstruct(G) 
    G = gsp_nn_graph(X, param.paramnn);
end

if ~isfield(G,'lmax')
    G = gsp_estimate_lmax(G);
end

%% Optimization
paraml1.verbose = param.verbose -1;
paraml1.y = X;
f1.prox = @(x,T) prox_l1(x,T,paraml1);
f1.eval = @(x) norm(x(:),1);

f2.grad = @(x) gamma*2*G.L*x;
f2.eval = @(x) gamma*sum(gsp_norm_tik(G,x));
f2.beta = 2*gamma*G.lmax;

Lr = solvep(X,{f1,f2},param);
Sp =X-Lr;

%% Optional output parameters
if nargout>3
    [U, S , V] = svdecon(Lr);
end


end


