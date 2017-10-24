function [finer_approximation]=gsp_pyramid_synthesis_single_interpolation_old(coarse_approximation,prediction_error,G,keep_inds,param)
%GSP_PYRAMID_SYNTHESIS_SINGLE Sythesize a single level of the graph pyramid transform 
%   Usage:  [next_coarse_approximation]=gsp_pyramid_synthesis_single(current_coarse_approximation,prediction_error,G,keep_inds,param);
%
%   Input parameters:
%         current_coarse_approximation : Coarse approximation of the signal on a reduced graph.
%         prediction_error             : Prediction error that was made when forming the current coarse approximation.
%         G                            : Graph structure on which the signal resides.
%         keep_inds                    : The indices of the vertices to keep when downsampling the graph and signal.
%         param                        : Contains optional additional parameters.
%   Output parameters:
%         next_coarse_approximation    : Coarse approximation of the signal on a higher resolution graph.
%   Additional parameters:
%         param.least_squares          : Set to 1 to use the least squares synthesis method.
%         param.use_landweber          : Set to 1 to use the Landweber iteration approximation in the least squares synthesis.
%         param.regularize_epsilon     : Interpolation parameter.
%         param.landweber_its          : Number of iterations in the Landweber approximation for least squares synthesis.
%         param.landweber_tau          : Parameter for the Landweber iteration.
%
%   'gsp_pyramid_synthesis_single(current_coarse_approximation,prediction_error,G,keep_inds,param)' 
%   synthesizes a single level of the graph pyramid transform.
%
%   See also:  
%
%   Demos:  
% 
%   References: 

%   AUTHOR : David I Shuman.
%   TESTING: 
%   REFERENCE:

if ~isfield(param, 'least_squares')
    least_squares=0;
else
    least_squares = param.least_squares;
end

if ~isfield(G,'lmax')
    G=gsp_estimate_lmax(G);
end

% Compute the next coarse approximation at a higher resolution level
if ~least_squares % i.e., use the direct synthesis method
    upsampled_coarse_approximation=zeros(G.N,1);
    upsampled_coarse_approximation(keep_inds)=coarse_approximation;
    prediction=gsp_interpolate_old(upsampled_coarse_approximation,G,keep_inds,param);
    finer_approximation=prediction+prediction_error;
else
    if ~isfield(param,'use_landweber')
        if (G.N>3000 || ~isfield(G,'e') || ~isfield(G,'U') )
            use_landweber=1;
        else
            use_landweber=0;
        end
    else
        use_landweber=param.use_landweber;
    end
    Nbar=length(keep_inds);
    if ~isfield(param, 'regularize_epsilon')
        regularize_epsilon=.005;
    else
        regularize_epsilon = param.regularize_epsilon;
    end
    if use_landweber
        if ~isfield(param, 'landweber_its')
            landweber_its=50;
        else
            landweber_its = param.landweber_its;
        end
        if ~isfield(param, 'landweber_tau')
            landweber_tau=1;
        else
            landweber_tau = param.landweber_tau;
        end
        % This is the Landweber iteration to apply the pseudoinverse 
        % (the efficiency of this set of operations can probably be improved)
        x_old=zeros(G.N,1);
        z=[coarse_approximation;prediction_error];
        PhiV1t=zeros(Nbar,G.N);
        green_kernel=@(x) 1./(x+regularize_epsilon);
        for i=1:Nbar
            PhiV1t(i,:)=(gsp_filter(G,green_kernel,gsp_delta(G.N,keep_inds(i)),param))';
        end
        for k=1:landweber_its
            [x_bar,y_bar]=gsp_pyramid_analysis_single_interpolation(x_old,G,keep_inds,param.h_filter,param);
            z_prime=[x_bar;y_bar];
            z_delt=z-z_prime;
            alpha_new=PhiV1t*z_delt((Nbar+1):end);
            x_up=zeros(G.N,1);
            x_up(keep_inds)=z_delt(1:Nbar);
            regularized_L= G.L+regularize_epsilon*eye(G.N);
            elim_inds=setdiff(1:G.N,keep_inds);
            next_term=regularized_L(keep_inds,keep_inds)*alpha_new-regularized_L(keep_inds,elim_inds)*(regularized_L(elim_inds,elim_inds)\(regularized_L(elim_inds,keep_inds)*alpha_new));
            next_up=zeros(G.N,1);
            next_up(keep_inds)=next_term;
            x_new=x_old+landweber_tau*(gsp_filter(G,param.h_filter,x_up-next_up,param)+z_delt((Nbar+1):end));
            x_old=x_new;
        end
        finer_approximation=x_new;
    else % When the graph is small enough, we can do a full eigendecomposition and compute the full analysis operator T_a
        if ( ~isfield(G,'e') || ~isfield(G,'U') )
            G=gsp_compute_fourier_basis(G);
        end
        H=G.U*diag(param.h_filter(G.e))*G.U';
        Phi=G.U*diag(1./(regularize_epsilon+G.e))*G.U'; 
        S=sparse(1:Nbar,keep_inds,ones(Nbar,1),Nbar,G.N);
        Ta=[S*H;eye(G.N)-Phi(:,keep_inds)*(Phi(keep_inds,keep_inds)\(S*H))];
        finer_approximation=(Ta'*Ta)\(Ta'*[coarse_approximation;prediction_error]);
    end
end
    
end