function [] = gsp_plot_jtv_filter(G, filters, filtertype, param)
%GSP_PLOT_JTV_FILTER  Plot a time-vertex filterbank
%   Usage:  gsp_plot_jtv_filter(G,filters);
%           gsp_plot_jtv_filter(G,filters,param);
%
%   Input parameters:
%       G          : Time-Vertex graph structure
%       filters    : Cell array of time-vertex filters
%       filtertype : Filter domain (ts,js,ts-array,js-array)
%       param      : Structure of optional parameters
%   Output parameters:
%       none
%
%   Example:
%
%         alpha = [0.1 0.5 1 2];
%         G = gsp_sensor(100);
%         G = gsp_jtv_graph(G,100,1);
%         G = gsp_estimate_lmax(G);
%         [g, filtertype] = gsp_jtv_design_wave(G, alpha);
%         param.domain='time-spectral';
%         gsp_plot_jtv_filter(G, g, filtertype,param);
%         param.domain='joint-spectral';
%         gsp_plot_jtv_filter(G, g, filtertype,param);
%
%
%   Additional parameters
%   ---------------------
%
%    param.npoints  : Number of points where the filters are evaluated if eigenvalues not available (default 100).
%    param.show_sum : Extra plot showing the sum of the squared magnitudes of the filters (default 1 if there is multiple filters).
%    param.verbose  : Verbosity level (1 display the warning - 0 no log) (default 1).
%    param.title    : Cell array of title for subplots. (default 1:Nf)
%    param.domain   : Visualize the spectrum in 'time-spectral' or 'joint-spectral' (default param.domain='joint-spectral')
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/plotting/gsp_plot_jtv_filter.php

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

% Author: Francesco Grassi, Nathanael Perraudin
% Date   : September 2016

% Read input parameters
if nargin < 4
    param = struct;
end

if nargin<3
    error('Invalid number of arguments: GSP_PLOT_JTV_FILTER needs the type of time-vertex filter.')
end

Nf=numel(filters);

if ~isfield(param,'show_sum'), param.show_sum = Nf>1; end
if ~isfield(param,'npoints'),  param.npoints = 100; end
if ~isfield(param,'title'),    param.title=num2cell(1:Nf); end
if ~isfield(param,'domain'),   param.domain='joint-spectral'; end
if ~isfield(param,'fftshift'), param.fftshift = 1;end
%% Define axis

if isfield(G,'e')
    lambdas = G.e;
else
    lambdas = linspace(0,G.lmax,param.npoints);
end

switch filtertype
    case  {'ts','ts-array'}
        v = gsp_jtv_ta(G);
    case  {'js','js-array'}
        v = gsp_jtv_fa(G);
    otherwise
        error('Unknown filtertype');
end


%% Evaluating filter

fid = gsp_jtv_filter_evaluate(filters,filtertype,lambdas,v,param);


if param.show_sum
    Nf = Nf+1;
    test_sum = sum(fid.^2,3);
    fid(:,:,Nf) = test_sum;
    param.title{Nf} = 'Sum squared coeff';
end


%% Plot

switch param.domain
    case 'time-spectral'
        [v,xlab] = gsp_jtv_ta(G);
        fid = real(fid);
    case 'joint-spectral'
        [v,xlab] = gsp_jtv_fa(G,param.fftshift);
        fid = abs(fid);
        if param.fftshift
            fid = fftshift(fid,2);
        end
    otherwise
        error('Unknown domain');
end



figure
r = max(1,floor(sqrt(Nf)));
c = max(1,ceil(Nf/r));
for ii=1:Nf
    subplot(r,c,ii)
    imagesc(v,lambdas,fid(:,:,ii))
    xlabel(xlab);
    ylabel('\lambda');
    title(param.title{ii})
    axis xy
end

end

