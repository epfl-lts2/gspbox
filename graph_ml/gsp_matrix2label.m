function [ f ] = gsp_matrix2label( B, minf )
%GSP_MATRIY2LABEL Reconstruct labels from matrix
%   Usage:  f = gsp_matrix2label( B );
%           f = gsp_matrix2label( B, minf );
%   
%   Input parameters:
%       B   : Classification matrix
%       minf: smallest integer (default 0)
%   Output parameters:
%       f   : Labels
%
%   See also: gsp_classification_matrix gsp_classification_tik gsp_classification_tv
%
%   Url: http://lts2research.epfl.ch/gsp/doc/graph_ml/gsp_matrix2label.php

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.5.1
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
% Date  : 24 July 2015

if nargin<2
    minf = 0;
end

[~,f] = max(B,[],2);
f = f-1+minf;

end


