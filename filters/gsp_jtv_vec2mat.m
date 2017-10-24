function X = gsp_jtv_vec2mat(G,X)
%GSP_JTV_VEC2MAT vector to matrix transform
%   Usage: d  = gsp_jtv_vec2mat(G,X);
%
%   Input parameters:
%       G       : Graph
%       X       : Data
%
%   Ouput parameter
%       X       : Data
%   
%   Reshape the data from the vector form to the matrix form


[Nx,Ny] = size(X);

X = reshape(X,G.N,G.jtv.T,Nx/G.N/G.jtv.T,Ny);


end