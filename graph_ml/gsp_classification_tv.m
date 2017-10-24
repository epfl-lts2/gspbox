function sol = gsp_classification_tv(G ,M, y , tau, param )
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

soltv = gsp_regression_tv(G, M, B , tau, param );

minf = min(y);
sol = gsp_matrix2label(soltv,minf);



end