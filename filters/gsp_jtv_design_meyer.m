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