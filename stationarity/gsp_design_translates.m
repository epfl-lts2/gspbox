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
%   This function construct a filter bank of *N* uniformly translated
%   filter from the mother filter *g0*.
%

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

