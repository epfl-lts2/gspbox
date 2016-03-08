
clear
close all;
[ x,y,xx,yy,f,ff ] = prepare_usps_full( );
%%[x, y, xx, yy] = load_usps_full();
%
%   Url: http://lts2research.epfl.ch/gsp/doc/graph_ml/classify_usps.php

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

% Nel = 300;
% x = x(:,1:Nel);
% xx = xx(:,1:Nel);
% y = y(1:Nel);
% yy = yy(1:Nel);
% f = f(:,1:Nel);
% ff = ff(:,1:Nel);
param.verbose = 1;
param.k = 6;

%%
%[rel_errorx, name ] = graph_mlcl_compare_all(x(:,1:1000)',xx(:,1:1000)',y(1:1000),yy(1:1000), param);
[rel_errorx, name ] = graph_mlcl_compare_all(x',xx',y,yy, param);
%%
%[rel_errorf, name ] = graph_mlcl_compare_all(f',ff',y,yy, param);

