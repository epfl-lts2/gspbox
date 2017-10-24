function [ g , mu ] = gsp_design_translates(G, g0,N )
%GSP_DESIGN_TRANSLATES Create a filterbank by uniformly translating a window
%   Usage: g = gsp_design_translates( G, g0, Ntrans );
%   
%   Inputs parameters:
%       G       : Graph structure
%       g0      : Mother window (anonymous function)
%       N       : Number of translate
%
%   Outputs parameters:
%       g       : filterbank
%       mu      : Centers of the filters
%
%   This function construct a filter bank of N uniformly translated
%   filter from the mother filter g0.
%
%
%   Url: https://epfl-lts2.github.io/gspbox-html/doc/stationarity/gsp_design_translates.html

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.7.4
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

% Author : Nathanael Perraudin
% Date: 6 January 2016


if isstruct(G)
    if ~isfield(G,'lmax')
        if param.verbose
            fprintf('GSP_DESIGN_TRANSLATE has to compute lmax \n')
        end
        G = gsp_estimate_lmax(G);
    end
   lmax = G.lmax;
else
   lmax = G;
end


mu = linspace(0,lmax,N);

g = cell(length(mu),1);

for ii = 1:length(mu)
    g{ii} = @(x) g0(x-mu(ii));
end


end


