% GSPBOX - Filters
%
%  Design - N filter
%    gsp_design_mexican_hat       -  Design a mexican hat filterbank
%    gsp_design_abspline          -  Design a abspline filterbank
%    gsp_design_meyer             -  Design a Meyer filterbank (tight)
%    gsp_design_simple_tf         -  Design a simple tight frame filterbank
%    gsp_design_itersine          -  Design a itersine filterbank (tight)
%    gsp_design_half_cosine       -  Design a half cosine filterbank (tight)
%    gsp_design_warped_translates -  Design a filterbank with a warping function 
%
%  Design - 2 filter (LP - HP) tight filterbank
%    gsp_design_regular           -  Design 2 filter with the "regular" construction
%    gsp_design_held              -  Design 2 filter with the "Held" construction
%    gsp_design_simoncelli        -  Design 2 filter with the "Simoncelli" construction
%    gsp_design_papadakis         -  Design 2 filter with the "Papadakis" construction
%
%  Application
%    gsp_filter                   -  Shortcut to gsp_filter_analysis
%    gsp_filter_analysis          -  Analysis operator for filterbank
%    gsp_filter_synthesis         -  Synthesis operator for filterbank
%    gsp_filter_inverse           -  Inverse operator for filterbank
%
%   Size Handling
%    gsp_mat2vec                  -  Matrix to vector representation for filterbanks 
%    gsp_vec2mat                  -  Vector to matrix representation for filterbanks 
%
%  Utils
%    gsp_filterbank_matrix        -  Create the analysis operator matrix associated to a filterbank 
%    gsp_wlog_scales              -  Compute log scale vector for wavelets
%    gsp_filter_evaluate          -  Evaluate a filterbank
%    gsp_filterbank_bounds        -  Bound for the filterbank
%    gsp_tighten_filter           -  Create a filter that tighten the filterbank
%
%  For help, bug reports, suggestions etc. please send email to
%  gspbox 'dash' support 'at' groupes 'dot' epfl 'dot' ch
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/filters/index.php

% Copyright (C) 2013-2014 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.3.1
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

% To be done
%   - add function in doc from meyer and simple_tf

