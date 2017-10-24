function [ d ] = gsp_vec2mat( d,Nf )
%GSP_VEC2MAT vector to matrix transform
%   Usage: d  = gsp_vec2mat( d,Nf );
%
%   Input parameters:
%       d       : Data
%       Nf      : Number of filter
%
%   Ouput parameter
%       d       : Data
%   
%   Reshape the data from the vector form to the matrix form

% TESTING: test_filter

[M,N] = size(d);

d = reshape(d,M/Nf,Nf,N);


end

