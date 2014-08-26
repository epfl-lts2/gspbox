function [c] = gsp_filter_analysis(G, fi, s, param)
%GSP_FILTER_ANALYSIS Analysis operator of a gsp filterbank
%   Usage:  coeffs = gsp_filter_analysis(G, fi, signal);
%           coeffs = gsp_filter_analysis(G, fi, signal, param);
%
%   Input parameters:
%         G         : Graph structure.
%         fi         : Set of spectral graph filters.
%         s         : graph signal to analyze.
%         param     : Optional parameter
%   Output parameters:
%         c         : Transform coefficients
%
%   'gsp_filter_analysis(G,fi,signal)' computes the transform
%   coefficients of a signal f, where the atoms of the transform
%   dictionary are generalized translations of each graph spectral filter
%   to each vertex on the graph.
%
%      c = D' * f 
%
%   where the columns of D are g_{i,m}=T_i g_m, and T_i is a
%   generalized translation operator applied to each filter 
%   hat{g}_m(cdot).  
%
%   Each column of c is the response of the signal to one filter.
%
%   Example:
%
%         Nf = 5;
%         param.distribute = 1;
%         G = gsp_sensor(256);
%         G = gsp_compute_fourier_basis(G);
%         paramf.log = 1;
%         g = gsp_design_warped_translates(G, Nf,paramf);  
%         s = sign(G.U(:,2));
%         sf = gsp_vec2mat(gsp_filter_analysis(G,g,s),Nf);
%         paramplot.show_edges = 1;
%         figure()
%         subplot(221)
%         gsp_plot_signal(G,sf(:,2),paramplot);
%         subplot(222)
%         gsp_plot_signal(G,sf(:,3),paramplot);       
%         subplot(223)
%         gsp_plot_signal(G,sf(:,4),paramplot);      
%         subplot(224)
%         gsp_plot_signal(G,sf(:,5),paramplot);
%       
%
%   Additional parameters
%   ---------------------
%  
%    param.exact  : To use exact graph spectral filtering instead of the
%     Chebyshev approximation. To use this option, the graph need the
%     Fourier basis of the graph need to be computed (default 0).
%    param.cheb_order : Degree of the Chebyshev approximation
%     (default=30). 
%    param.verbose : Verbosity level (0 no log - 1 display warnings)
%     (default 1).   
%
%   See also: gsp_filter_synthesis gsp_filter_inverse
% 
%   References:
%     D. K. Hammond, P. Vandergheynst, and R. Gribonval. Wavelets on graphs
%     via spectral graph theory. Appl. Comput. Harmon. Anal., 30(2):129-150,
%     Mar. 2011.
%     
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/filters/gsp_filter_analysis.php

% Copyright (C) 2013-2014 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.3.0
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

% Author: David I Shuman, Nathanael Perraudin
% Testing: test_filter
% Date: 19 March 2014

  
% Read input parameters
if nargin < 4
    param = struct;
end



Nf = length(fi);

if ~isfield(param,'exact'); param.exact = isfield(G,'U'); end
if ~isfield(param,'cheb_order'); param.cheb_order = 30; end
if ~isfield(param,'verbose'); param.verbose = 1; end



if param.exact 
    if ( ~isfield(G,'e') || ~isfield(G,'U') )
        if param.verbose
            warning(['GSP_FILTER_ANALYSIS: The Fourier matrix is not ',...
                'available. The function will compute it for you. ',...
                'However, if you apply many time this function, you ',...
                'should precompute it using the function: ',...
                'gsp_compute_fourier_basis']);
        end
        G=gsp_compute_fourier_basis(G);
    end
    Nv = size(s,2);
    c = zeros(G.N*Nf,Nv);
    
    fie = gsp_filter_evaluate(fi,G.e);
    
    for ii=1:Nf
        c((1:G.N)+G.N * (ii-1),:)= gsp_igft(G, ...
            repmat(fie(:,ii),1,Nv) ...
            .* gsp_gft(G, s));
    end


else % use Chebyshev approximation method


    
    if ~isfield(G,'lmax');
        G = gsp_estimate_lmax(G);
        if param.verbose
            warning(['GSP_FILTER_ANALYSIS: The variable lmax is not ',...
                'available. The function will compute it for you. ',...
                'However, if you apply many time this function, you ',...
                'should precompute it using the function: ',...
                'gsp_estimate_lmax']);
        end
    end
    

    cheb_coeffs = gsp_cheby_coeff(G, fi,...
            param.cheb_order, param.cheb_order +1);      
    c = gsp_cheby_op(G, cheb_coeffs, s);


end


