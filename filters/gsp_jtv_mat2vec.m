function C = gsp_jtv_mat2vec(C)
%GSP_JTV_MAT2VEC vector to matrix transform
%   Usage:  C  = gsp_jtv_mat2vec(C);
%
%   Input parameters:
%       C       : Data
%
%   Ouput parameter
%       C       : Data
%   
%   Reshape the data from the matrix form to the vector form

[N,T,Nf,Ns] = size(C);


C = reshape(C,N*T*Nf,Ns);

end