function [finer_approximation]=gsp_pyramid_synthesis_single(G,coarse_approximation,prediction_error,keep_inds,param)
%GSP_PYRAMID_SYNTHESIS_SINGLE Synthesize a single level of the graph pyramid transform 
%   Usage:  [finer_approximation]=gsp_pyramid_synthesis_single(G,coarse_approximation,prediction_error,keep_inds,param);
%
%   Input parameters:
%         G                    : Graph structure on which the signal resides.
%         coarse_approximation : Coarse approximation on a reduced graph.
%         prediction_error     : Prediction error that was made when forming the coarse approximation.
%         keep_inds            : Indices of the vertices to keep when downsampling the graph and signal.
%         param                : Additional parameters.
%   Output parameters:
%         finer_approximation  : Coarse approximation of the signal on a higher resolution graph.
%
%   'gsp_pyramid_synthesis_single(G,coarse_approximation,prediction_error,keep_inds,param)' 
%   synthesizes a single level of the graph pyramid transform.
%  
%   param is a structure containing optional arguments including
%
%    param.regularize_epsilon : Interpolation parameter.
%    param.least_squares : Set to 1 to use the least squares synthesis 
%     (default=0)  
%    param.use_landweber : Set to 1 to use the Landweber iteration 
%     approximation in the least squares synthesis.
%    param.landweber_its : Number of iterations in the Landweber 
%     approximation for least squares synthesis.
%    param.landweber_tau : Parameter for the Landweber iteration.
%    param.h_filter : A graph spectral filter. This filter is 
%     required for least squares synthesis, but not for the direct 
%     synthesis method 
%
%   Please read the documentation of GSP_FILTER_ANALYSIS for other
%   optional arguments.
%
%   See also: gsp_graph_multiresolution gsp_pyramid_synthesis 
%             gsp_pyramid_analysis gsp_pyramid_analysis_single 
%             gsp_pyramid_cell2coeff
%
%   Demo: gsp_demo_pyramid
% 
%   References:
%     D. I. Shuman, M. J. Faraji, and P. Vandergheynst. A framework for
%     multiscale transforms on graphs. arXiv preprint arXiv:1308.4942, 2013.
%     
%     
%
%   Url: http://lts2research.epfl.ch/gsp/doc/operators/gsp_pyramid_synthesis_single.php

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.6.0
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

%   Author : David I Shuman, Nathanael Perraudin.
%   Date   : 26 November 2015
%   Testing: test_pyramid

if ~isfield(param, 'least_squares')
    param.least_squares=0;
end

if ~isfield(G,'lmax')
    G=gsp_estimate_lmax(G);
end

% Compute the next coarse approximation at a higher resolution level
if ~param.least_squares % i.e., use the direct synthesis method
%     upsampled_coarse_approximation=zeros(G.N,size());
%     upsampled_coarse_approximation(keep_inds,:)=coarse_approximation;
    
    prediction=gsp_interpolate(G,coarse_approximation,keep_inds,param);
    finer_approximation=prediction+prediction_error;
else
    if ~isfield(param, 'h_filter')
        error('h-filter not provided');
    end
    if ~isfield(param,'use_landweber')
        if (G.N>3000 || ~isfield(G,'e') || ~isfield(G,'U') )
            param.use_landweber=1;
        else
            param.use_landweber=0;
        end
    end
    Nbar=length(keep_inds);
    if ~isfield(param, 'regularize_epsilon')
        param.regularize_epsilon=.005;
    end
    if param.use_landweber
        if ~isfield(param, 'landweber_its')
            param.landweber_its=50;
        end
        if ~isfield(param, 'landweber_tau')
            param.landweber_tau=1;
        end
        % This is the Landweber iteration to apply the pseudoinverse 
        % (the efficiency of this set of operations can probably be improved)
        x_old=zeros(G.N,1);
        z=[coarse_approximation;prediction_error];
        PhiV1t=zeros(Nbar,G.N);
        green_kernel=@(x) 1./(x+param.regularize_epsilon);
        for i=1:Nbar
            PhiV1t(i,:)=(gsp_filter(G,green_kernel,gsp_delta(G.N,keep_inds(i)),param))';
        end
        for k=1:param.landweber_its
            [x_bar,y_bar]=gsp_pyramid_analysis_single(G,x_old,keep_inds,param.h_filter,param);
            z_prime=[x_bar;y_bar];
            z_delt=z-z_prime;
            alpha_new=PhiV1t*z_delt((Nbar+1):end);
            x_up=zeros(G.N,1);
            x_up(keep_inds)=z_delt(1:Nbar);
            regularized_L= G.L+param.regularize_epsilon*eye(G.N);
            elim_inds=setdiff(1:G.N,keep_inds);
            next_term=regularized_L(keep_inds,keep_inds)*alpha_new-regularized_L(keep_inds,elim_inds)*(regularized_L(elim_inds,elim_inds)\(regularized_L(elim_inds,keep_inds)*alpha_new));
            next_up=zeros(G.N,1);
            next_up(keep_inds)=next_term;
            x_new=x_old+param.landweber_tau*(gsp_filter(G,param.h_filter,x_up-next_up,param)+z_delt((Nbar+1):end));
            x_old=x_new;
        end
        finer_approximation=x_new;
    else % When the graph is small enough, we can do a full eigendecomposition and compute the full analysis operator T_a
        if ( ~isfield(G,'e') || ~isfield(G,'U') )
            G=gsp_compute_fourier_basis(G);
        end
        H=G.U*diag(param.h_filter(G.e))*G.U';
        Phi=G.U*diag(1./(param.regularize_epsilon+G.e))*G.U'; 
        S=sparse(1:Nbar,keep_inds,ones(Nbar,1),Nbar,G.N);
        Ta=[S*H;eye(G.N)-Phi(:,keep_inds)*(Phi(keep_inds,keep_inds)\(S*H))];
        finer_approximation=(Ta'*Ta)\(Ta'*[coarse_approximation;prediction_error]);
    end
end
    
end
