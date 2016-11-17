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
%
%   Url: http://lts2research.epfl.ch/gsp/doc/utils/gsp_repmatline.php

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.7.0
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


