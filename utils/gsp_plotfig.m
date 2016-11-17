function [ ] = gsp_plotfig( save_name,param )
%GSP_PLOTFIG Ploting figures with optimal size for paper
%   Usage: gsp_plotfig(save_name);
%          gsp_plotfig(save_name,param);
%   
%   Input parameters:
%       save_name: name to save the figure
%       param   : optional parameters
%
%   param a Matlab structure containing the following fields:
%
%    param.pathfigure : path to the folder to save the figures
%    param.legendlocation : location of the figure (default 'Best');
%    param.position : position and size of the figure 
%     (default [100 100 600 400])
%    param.labelsize : Size of the label (default 12)
%    param.titlesize : Size of the title (default 16)
%    param.titleweight*: Weight of the title (default 'normal')
%    param.save*: Save the figure (default 1)
%    param.eps*: Save the figure in eps instead of png (default 0)
%
%   Url: http://lts2research.epfl.ch/gsp/doc/utils/gsp_plotfig.php

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


% Nathanael Perraudin
% 19 June 2013

if nargin<2
    param=struct;
end

% Optional parameters
if ~isfield(param, 'pathfigure'), param.pathfigure = 'figures/'; end
if ~isfield(param, 'position'), param.position = [100 100 600 400]; end
if ~isfield(param, 'labelsize'), param.labelsize = 14; end
if ~isfield(param, 'titlesize'), param.titlesize = 16; end
if ~isfield(param, 'titleweight'), param.titleweight = 'normal'; end
if ~isfield(param, 'save'), param.save =1; end
% if ~isfield(param, 'baw'), param.baw =0; end
if ~isfield(param, 'eps'), param.eps =0; end

    


% set the axes
try  %#ok<TRYNC>
    set(gcf, 'Position', param.position);
    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperSize',[param.position(3)*1.02, param.position(4)*1.02]);
end


allAxes = findall(0,'type','axes');
for ii = 1:numel(allAxes)
    % set the title
    try %#ok<TRYNC>
        set(allAxes(ii),'FontSize',param.labelsize)
    end
    try %#ok<TRYNC>
        h=get(allAxes(ii),'Title');
        set(h,'FontSize',param.titlesize);
        set(h,'FontWeight',param.titleweight);
    end
    try %#ok<TRYNC>
        h=get(allAxes(ii),'xlabel');
        set(h,'FontSize',param.labelsize);
    end
    try %#ok<TRYNC>
        h=get(allAxes(ii),'ylabel');
        set(h,'FontSize',param.labelsize);
    end
end
drawnow;

% save the results



if param.save
    if ~isdir(param.pathfigure)
           mkdir(param.pathfigure);
    end
    filename=strcat(param.pathfigure,save_name);
    if param.eps
         %print('-depsc',[filename,'.eps']);
         print2eps([filename,'.eps']);
    else
        %print('-dpng','-opengl','-r300',[filename,'.png']);
        print('-dpng','-r300',[filename,'.png']);
        %hgsave([filename,'.fig']);
    end
end
end


