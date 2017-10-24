function [ Ar ] = gsp_repmatline( A,ncol,nrow )
%GSP_REPMATLINE This function repeat the matrix A in a specific manner
%   Usage Ar = gsp_repmatline( A,ncol,nrow );
%
%   Inputs parameters
%       A   : Matrix
%       ncol: Integer
%       nrow: Integer
%   Outputs parameters
%       Ar  : Matrix
%
%   This function repeat a matrix line by line and column by column
%
%   For ncol=1 and nrow=2, the matix 
%               1 2
%               3 4
%   becomes
%               1 1 2 2
%               3 3 4 4
%

% Author: Nathanael Perraudin
% Date  : 23 October 2013

if ~floor(ncol)==ncol || ~floor(nrow)==nrow
   error('The number of lines and rows must be integer');
end

if ncol<1 || nrow<1
   error('The number of lines and rows must be greater or equal to one') 
end


% Initialisation
Ar=A;

if ncol>1

    [row, col] = size(Ar);
    Ar = Ar(:);
    Ar = Ar(:,ones(1,ncol)).';
    Ar = reshape(Ar,ncol*row,col);
end

if nrow>1
    Ar=transpose(Ar);
    [col, row] = size(Ar);
    Ar = Ar(:);
    Ar = Ar(:,ones(1,nrow)).';
    Ar = reshape(Ar,nrow*col,row);
    Ar=transpose(Ar);
end


end

