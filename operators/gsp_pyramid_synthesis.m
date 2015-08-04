function [signal,ca] = gsp_pyramid_synthesis(Gs,coeff,param)
%GSP_PYRAMID_SYNTHESIS Synthesizes a signal from its graph pyramid transform coefficients 
%   Usage:  signal = gsp_pyramid_synthesis(Gs,coeff)
%           [signal, ca ] = gsp_pyramid_synthesis(Gs,coeff)
%
%   Input parameters:
%       Gs      : A multiresolution sequence of graph structures.
%       coeff   : The coefficients to perform the reconstruction
%   Output parameters:
%       signal  : The synthesized signal.
%       ca      : Cell array with the coarse approximation at each level
%
%   This function perform the pyramid synthesis of the coefficient in
%   coeff. 
%
%   The pyramid analysis operator returns two arguments:
%
%           [ca,pe]=gsp_pyramid_analysis(Gs, f);
%
%   To obtain the coefficients you can call the function:
%
%           coeff = gsp_pyramid_cell2coeff(ca,pe);
%
%   Example:
%
%       N = 256;
%       Nl = 4;
%       G = gsp_sensor(N);
%       Gs = gsp_kron_pyramid(G,Nl);
%       Gs = gsp_compute_fourier_basis(Gs);
%       f = rand(N,1);
%       [ca,pe]=gsp_pyramid_analysis(Gs, f);
%       coeff = gsp_pyramid_cell2coeff(ca,pe);
%       f_pred = gsp_pyramid_synthesis(Gs,coeff);
%       error = norm(f-f_pred)
%
%   See also: gsp_pyramid_analysis gsp_kron_pyramid gsp_pyramid_cell2coeff
%
%   Demo: gsp_demo_pyramid
%
%   References:
%     I. Pesenson. Variational splines and paley--wiener spaces on
%     combinatorial graphs. Constructive Approximation, 29(1):1--21, 2009.
%     
%     D. I. Shuman, M. J. Faraji, and P. Vandergheynst. A framework for
%     multiscale transforms on graphs. arXiv preprint arXiv:1308.4942, 2013.
%     
%
%   Url: http://lts2research.epfl.ch/gsp/doc/operators/gsp_pyramid_synthesis.php

% Copyright (C) 2013-2014 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.4.0
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

% Author: Nathanael Perraudin
% Date  : 5 August 2014
% Testing : test_pyramid


if nargin < 3
    param = struct;
end
    
if ~isfield(param,'order'), param.order = 100; end

Nl = length(Gs)-1;

% Initisalization
Nt = Gs{Nl+1}.N;
ca{Nl+1} = coeff(1:Nt, :);

ind = Nt+1;
% Reconstruct each level
for ii = 1:Nl
      Nt = Gs{Nl+1-ii}.N;
    % Compute prediction
    s_pred = gsp_interpolate(Gs{Nl+1-ii},Gs{Nl+2-ii},ca{Nl+2-ii},param);
    % Compute the ca coeff
    ca{Nl+1-ii} = s_pred + coeff(ind:(ind+Nt-1), :);
    
    ind = ind + Nt;
end

signal = ca{1};

end





% Here below is David code that might be used if we want to implement
% more... Right now. we just do it simple


% 
% 
% % Read input parameters 
% if nargin>4
%     param=varargin{1};
% else
%     param=0;
% end
% 
% num_levels=length(prediction_errors);
% 
% % Compute the pyramid transform
% current_coarse_approximation=coarsest_approximation;
% 
% for i=num_levels:-1:1
%     if isfield(param, 'least_squares')
%         if param.least_squares
%             if ~isfield(param, 'h_filters')
%                 error('h-filter not provided');
%             else
%                 param.h_filter=param.h_filters{i};
%             end
%         end
%     end
%     current_coarse_approximation=gsp_pyramid_synthesis_single_interpolation(current_coarse_approximation,prediction_errors{i},Gs{i},subsampled_vertex_indices{i},param);
% end
% reconstruction=current_coarse_approximation;
% 
% function [finer_approximation]=gsp_pyramid_synthesis_single_interpolation(coarse_approximation,prediction_error,G,keep_inds,param)
% %GSP_PYRAMID_SYNTHESIS_SINGLE Sythesize a single level of the graph pyramid transform 
% %   Usage:  [next_coarse_approximation]=gsp_pyramid_synthesis_single(current_coarse_approximation,prediction_error,G,keep_inds,param);
% %
% %   Input parameters:
% %         current_coarse_approximation : Coarse approximation of the signal on a reduced graph.
% %         prediction_error             : Prediction error that was made when forming the current coarse approximation.
% %         G                            : Graph structure on which the signal resides.
% %         keep_inds                    : The indices of the vertices to keep when downsampling the graph and signal.
% %         param                        : Contains optional additional parameters.
% %   Output parameters:
% %         next_coarse_approximation    : Coarse approximation of the signal on a higher resolution graph.
% %   Additional parameters:
% %         param.least_squares          : Set to 1 to use the least squares synthesis method.
% %         param.use_landweber          : Set to 1 to use the Landweber iteration approximation in the least squares synthesis.
% %         param.regularize_epsilon     : Interpolation parameter.
% %         param.landweber_its          : Number of iterations in the Landweber approximation for least squares synthesis.
% %         param.landweber_tau          : Parameter for the Landweber iteration.
% %
% %   'gsp_pyramid_synthesis_single(current_coarse_approximation,prediction_error,G,keep_inds,param)' 
% %   synthesizes a single level of the graph pyramid transform.
% %
% %   See also:  
% %
% %   Demos:  
% % 
% %   References: 
% 
% %   AUTHOR : David I Shuman.
% %   TESTING: 
% %   REFERENCE:
% 
% if ~isfield(param, 'least_squares')
%     least_squares=0;
% else
%     least_squares = param.least_squares;
% end
% 
% % Compute the next coarse approximation at a higher resolution level
% if ~least_squares % i.e., use the direct synthesis method
%     upsampled_coarse_approximation=zeros(G.N,1);
%     upsampled_coarse_approximation(keep_inds)=coarse_approximation;
%     prediction=gsp_interpolate(upsampled_coarse_approximation,G,keep_inds,param);
%     finer_approximation=prediction+prediction_error;
% else
%     if ~isfield(param,'use_landweber')
%         if (G.N>3000 || ~isfield(G,'E') || ~isfield(G,'U') )
%             use_landweber=1;
%         else
%             use_landweber=0;
%         end
%     else
%         use_landweber=param.use_landweber;
%     end
%     Nbar=length(keep_inds);
%     if ~isfield(param, 'regularize_epsilon')
%         regularize_epsilon=.005;
%     else
%         regularize_epsilon = param.regularize_epsilon;
%     end
%     if use_landweber
%         if ~isfield(param, 'landweber_its')
%             landweber_its=50;
%         else
%             landweber_its = param.landweber_its;
%         end
%         if ~isfield(param, 'landweber_tau')
%             landweber_tau=1;
%         else
%             landweber_tau = param.landweber_tau;
%         end
%         % This is the Landweber iteration to apply the pseudoinverse 
%         % (the efficiency of this set of operations can probably be improved)
%         x_old=zeros(G.N,1);
%         z=[coarse_approximation;prediction_error];
%         PhiV1t=zeros(Nbar,G.N);
%         green_kernel=@(x) 1./(x+regularize_epsilon);
%         for i=1:Nbar
%             PhiV1t(i,:)=(gsp_filter(sgwt_delta(G.N,keep_inds(i)),G,green_kernel,param))';
%         end
%         for k=1:landweber_its
%             [x_bar,y_bar]=gsp_pyramid_analysis_single_interpolation(x_old,G,keep_inds,param.h_filter,param);
%             z_prime=[x_bar;y_bar];
%             z_delt=z-z_prime;
%             alpha_new=PhiV1t*z_delt((Nbar+1):end);
%             x_up=zeros(G.N,1);
%             x_up(keep_inds)=z_delt(1:Nbar);
%             regularized_L= G.L+regularize_epsilon*eye(G.N);
%             elim_inds=setdiff(1:G.N,keep_inds);
%             next_term=regularized_L(keep_inds,keep_inds)*alpha_new-regularized_L(keep_inds,elim_inds)*(regularized_L(elim_inds,elim_inds)\(regularized_L(elim_inds,keep_inds)*alpha_new));
%             next_up=zeros(G.N,1);
%             next_up(keep_inds)=next_term;
%             x_new=x_old+landweber_tau*(gsp_filter(x_up-next_up,G,param.h_filter,param)+z_delt((Nbar+1):end));
%             x_old=x_new;
%         end
%         finer_approximation=x_new;
%     else % When the graph is small enough, we can do a full eigendecomposition and compute the full analysis operator T_a
%         if ( ~isfield(G,'E') || ~isfield(G,'U') )
%             G=gsp_full_eigen(G);
%         end
%         H=G.U*diag(param.h_filter(G.E))*G.U';
%         Phi=G.U*diag(1./(regularize_epsilon+G.E))*G.U'; 
%         S=sparse(1:Nbar,keep_inds,ones(Nbar,1),Nbar,G.N);
%         Ta=[S*H;eye(G.N)-Phi(:,keep_inds)*pinv(Phi(keep_inds,keep_inds))*S*H];
%         finer_approximation=(Ta'*Ta)\(Ta'*[coarse_approximation;prediction_error]);
%     end
% end
%     
% end

