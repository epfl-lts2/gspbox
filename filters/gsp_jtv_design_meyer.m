function [g,filtertype] = gsp_jtv_design_meyer(G,Nf)
%GSP_DESIGN_JTV_MEYER Design the jtv Meyer tight filterbank
%   Usage: [g,filtertype] = gsp_design_jtv_meyer(G);
%          [g,filtertype] = gsp_design_jtv_meyer(G,Nf);
%
%   Input parameters:
%       G       : Time-Vertex graph structure
%       Nf      : Number of filters for each domain (total number Nf^2, default Nf = 4)
%   Output parameters:
%       g          : Cell array of time-vertex filters
%       filtertype : Filter domain js
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/filters/gsp_jtv_design_meyer.php

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

% Author :  Francesco Grassi
% Date   : September 2016

if nargin<2
    Nf = 4;
end

%graph meyer
g1 = gsp_design_meyer(G,Nf);

%time meyer
g2 = gsp_design_meyer(0.5/G.jtv.fs,Nf);

%building jtv meyer separable filterbank
n = 0;
g = cell(Nf^2,1);
for ii=1:Nf
    for jj=1:Nf
        n = n+1;
        g{n} = @(x,y) g1{ii}(x).*g2{jj}(abs(y));
    end
end

filtertype = 'js';

end
