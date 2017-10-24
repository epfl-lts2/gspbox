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
%   Lr by solving an optimization problem:
%
%      argmin_Lr || Lr - X ||_1  + gamma tr( X^T L X) 
%
%   The sparse part is given by  S = X - L_r . 
%   
%   If 0 is given for G, the corresponding graph will
%   be computed internally. The graph construction can be tuned using the
%   optional parameter: param.paramnn.
%
%   Other optional parameters are used in optimization. For details about
%   those, please see the help of the function SOLVEP . 
%
%   If the number of output argument is greater than 2. The function, will
%   additionally compute a very economical SVD such that  Lr = U S V^T.
%
%   This function uses the UNLocBoX to be working.
%
%   References:
%     N. Shahid, N. Perraudin, V. Kalofolias, and P. Vandergheynst. Fast
%     robust pca on graphs. arXiv preprint arXiv:1507.08173, 2015.
%     
%     
%
%   See also: gsp_frpcaog_2g
%
%   Url: https://epfl-lts2.github.io/gspbox-html/doc/graph_ml/gsp_frpcaog_1g.html

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.7.4
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% If you use this toolbox please kindly cite
%     N. Perraudin, J. Paratte, D. Shuman, V. Kalofolias, P. Vandergheynst,
%     and D. K. Hammond. GSPBOX: A toolbox for signal processing on graphs.
%     ArXiv e-prints, Aug. 2014.
% http://arxiv.org/abs/1408.5781

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



