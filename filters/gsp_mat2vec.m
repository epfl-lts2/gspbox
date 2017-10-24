function [ d ,Nf] = gsp_mat2vec( d )
%GSP_MAT2VEC vector to matrix transform
%   Usage:  d  = gsp_mat2vec( d );
%          [ d ,Nf] = gsp_mat2vec( d );
%
%   Input parameters:
%       d       : Data
%
%   Ouput parameter
%       d       : Data
%       Nf      : Number of filter
%   
%   Reshape the data from the matrix form to the vector form

% TESTING: test_filter

if iscell(d)
    Nc = numel(d);
    d2 = cell(Nc,1);
    for ii = 1:Nc
        d2{ii} = gsp_mat2vec(d{ii});
    end
    d = d2;
    return
end

[M,Nf,N] = size(d);

d = reshape(d,M*Nf,N);


end

