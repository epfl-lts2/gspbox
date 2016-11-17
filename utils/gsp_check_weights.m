function [ a ] = gsp_check_weights( W )
%GSP_CHECK_WEIGHTS Check a weight matrix
%   Usage: a = gsp_check_weights( W );
%
%   Input parameters:
%       W       : Weight matrix
%   Output parameters:
%       a       : Warning code
%   
%   This function performs various test on the weight matrix W. It
%   returns:
%
%        0     : Everything is ok
%        1     : The martrix contains inf values
%        2     : The diagonal is not  0
%        3     : The matrix is not square
%        4     : The matrix contains nan values
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/utils/gsp_check_weights.php

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
% Date  : 12 june 2014

a = 0;


if sum(sum(isinf(W)))
    warning(['GSP_CHECK_WEIGHTS: There is infinite value',...
             ' in the weight matrix']);
    a = 1;
end

if sum(abs(diag(W)))
    disp(['GSP_CHECK_WEIGHTS: The diagonal',...
             ' of the weight matrix is not 0!']);
    a = 2;
end

if size(W,1) ~= size(W,2)
    warning('GSP_CHECK_WEIGHTS: The weight matrix is not square!');
    a = 3;
end

if sum(sum(isnan(W)))
    warning(['GSP_CHECK_WEIGHTS: There is infinite value',...
             ' in the weight matrix']);
    a = 4;
end

end


