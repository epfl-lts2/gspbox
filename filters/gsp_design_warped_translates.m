function [ g, wf ] = gsp_design_warped_translates( G, Nf, param )
%GSP_DESIGN_WARPED_TRANSLATES Create a vertex frequency filterbank
%   Usage: g = gsp_design_warped_translates( G, Nf );
%          g = gsp_design_warped_translates( G, Nf, param );
%   
%   Inputs parameters:
%       G       : Graph structure (or lmax for 'log' and 'log_plus' only) 
%       Nf      : Number of filter 
%       param   : Structure of optional parameters
%
%   Outputs parameters:
%       g       : filterbanks
%       wf      : warped function
%
%   This function designs filters that are warped versions of the uniform
%   half cosine translates described above.
%
%   This function will compute the maximum eigenvalue of the laplacian. To
%   be more efficient, you can precompute it using::
%
%       G = gsp_estimate_lmax(G);
%
%   Example:::
%
%             figure();
%             Nf = 10;
%             G = gsp_sensor(100);
%             G = gsp_estimate_lmax(G);
%             G = gsp_spectrum_cdf_approx(G);
%             g = gsp_design_warped_translates(G, Nf);   
%             gsp_plot_filter(G,g);
%             [A,B] = gsp_filterbank_bounds(G,g)
%
%   *param* is an optional structure containing the following fields
%
%   * *param.verbose*: verbosity level. 0 no log - 1 display warnings.
%     (default 1) 
%   * *param.warping_type*: Create a warping function according two
%     different methods (default 'spectrum_approximation'). Please read
%     below for more information about this parameter.
%   * *param.log*: On top of the other warping add a log function. An
%     alternative way to construct spectral graph wavelets. These are
%     adapted to the specific spectrum, not just the length of 
%     the spectrum. The final warping function will be:
%     
%     ..    log(f(x))
%
%     .. math:: \log(f(x))
%
%     where the function $f(x)$ is defined by the attribute
%     *param.warping_type*.
%     Warning: Additional required inputs: *param.warp_function*.
%   * *param.warp_function*: To provide a special warping function. This
%     parameter is used when *param.warping_type* is 'custom'.
%   * *param.interpolation_type*: select the interpolation type for the
%     spectrum samples. You can choose 'pwl' (piece wise linear)
%     or 'monocubic'. This attribute is used only when *param.warping_type* 
%     is 'spectrum_interpolation'. (default 'monocubic')
%   * *param.filter*: select the initial uniform filterbank 'half_cosine'
%     or 'itersine'. See gsp_design_half_cosine and
%     gsp_design_itersine for more information about those filterbank.
%     If you want to use your personal filter, just put it there. For
%     instance::
%               
%               param.filter = gsp_design_abspline(G,Nf);
%
%     (Default 'half_cosine')
%
%   * *param.overlap*: overlap of the initial filter. Works only
%     with *param.filter* set to 'itersine'. For tight frame, input an even
%     number (default 2).
%   
%
%   Warping methods
%   ---------------
%
%   The different warping type available in *param.warping_type* are:
%
%   * 'spectrum_interpolation': Warping functions based on spectrum
%     samples. From the samples, an approximation of the spectrum cdf is
%     obtained by interpolation. Then this function is used for the
%     warping. (i.e., like the filter banks [1] in Section 2, these are
%     spectrum-adapted filter banks).
%     
%     If you use this method you need to specify the input
%     *param.approx_spectrum* that contains two fields:
%     *param.approx_spectrum.x* and *param.approx_spectrum.y* that are the
%     the point of the cumulative density distribution.
%
%   * 'spectrum_approximation': This function will compute an approximation
%     of the cumulative density function of the graph Laplacian
%     eigenvalues and use it as warping function. If you want to use the
%     cdf later, you should precompute it using::
%
%       G = gsp_spectrum_cdf_approx(G);
%   
%   * 'custom': The user provide the warping function in the parameter:
%     *param.warp_function*.
%
%   References: shuman2013spectrum

% Author: Nathanael Perraudin, David Shuman
% Date  : 18 June 2014
% Testing: test_filter

% TODO: check the log case!

if nargin < 3
    param = struct;
end

if ~isfield(param, 'verbose'), param.verbose = 1; end;
if ~isfield(param, 'warping_type')
    param.warping_type = 'spectrum_approximation'; 
end;
if ~isfield(param, 'log'), param.log = 0; end;
if ~isfield(param, 'filter'), param.filter = 'half_cosine'; end
if ~isfield(param, 'overlap'), param.overlap = 2; end
if ~isfield(param, 'logmax'), param.logmax = 10; end
if ~isfield(param, 'interpolation_type')
    param.interpolation_type = 'monocubic';
end





if isstruct(G)
    if ~isfield(G,'lmax')
        if param.verbose
            fprintf('GSP_DESIGN_WARPED_TRANSLATES: has to compute lmax \n')
        end
        G = gsp_estimate_lmax(G);
    end
   lmax = G.lmax;
else
   lmax = G;
end







switch param.warping_type


    case 'spectrum_interpolation'
        
        if ~isfield(param,'approx_spectrum')
        	error('GSP_DESIGN_WARPED_TRANSLATES: Approximate spectrum is not specified')
        end
    
        switch param.interpolation_type
        
            case 'pwl'
               % error('To be debugged!')
                if ~isfield(param,'warp_start_one') % DIS: Get rid of two cases - just include in interpolation points
                    warp_start_one=0;
                else
                    warp_start_one=param.warp_start_one;
                end

                % Generate uniform translates covering [0,upper_bound_translates] % DIS: Get rid of two cases - just include in interpolation points
                if warp_start_one
                    wf = @(s) gsp_pwl_warp_fn(param.approx_spectrum.x,param.approx_spectrum.y,s);
                   % upper_bound_translates=max(param.approx_spectrum.y);
                else
                    wf = @(s) gsp_pwl_warp_fn(param.approx_spectrum.x,param.approx_spectrum.y,s);
                   % upper_bound_translates=max(param.approx_spectrum.y)-1;
                end
            case 'monocubic'

                % Generate uniform translates covering [0,upper_bound_translates]
                wf = @(s) gsp_mono_cubic_warp_fn(param.approx_spectrum.x,param.approx_spectrum.y,s);
            otherwise
            
            error('GSP_DESIGN_WARPED_TRANSLATES: Unknown interpolation type')
        end
            
    
    

        
    case 'spectrum_approximation'
         if ~isfield(G,'spectrum_cdf_approx')
            if param.verbose
                fprintf(['GSP_DESIGN_WARPED_TRANSLATES: has to compute',...
                    ' the spectrum continuous density function ',...
                    'approximation \n']);
            end
            G = gsp_spectrum_cdf_approx(G);
         end
        if isfield(param,'warp_function')
            error('GSP_DESIGN_WARPED_TRANSLATES: Custom warp function is defined but not used')
        else
            wf=G.spectrum_cdf_approx;
        end

%         if ~isfield(param,'upper_bound_translates')
%             param.upper_bound_translates=wf(lmax);
%             %error('Upper bound of the uniform translates is not specified')
%         end

    %%%%%%%%%%%%%%%%%
    % Other user-specified warping function
    case 'custom'


        if ~isfield(param,'warp_function')
            fprintf('GSP_DESIGN_WARPED_TRANSLATES: Custom warp function is not specified (Taking a uniform warping)\n')
            wf = @(x) x;
        else
            wf = param.warp_function;
        end

%         if ~isfield(param,'upper_bound_translates')
%             param.upper_bound_translates=wf(lmax);
%             %error('Upper bound of the uniform translates is not specified')
%         end

    otherwise
        error('GSP_DESIGN_WARPED_TRANSLATES: Warp function type not recognized')
end



if param.log
    xmax = wf(lmax);
    wflog = @(s) log(1+wf(s)/xmax*lmax*param.logmax+eps);
    wf = wflog;
end

% Old way of doing the log...
%     % Generate (num_filters-1) uniform translates covering [0,log(lmax)] 
%     switch param.filter
%         case 'half_cosine'
%             uniform_filters = gsp_design_half_cosine(wflog(lmax),Nf-1);
%         case 'itersine'
%             paramt.overlap = param.overlap;
%             uniform_filters = gsp_design_itersine(wflog(lmax),Nf-1,paramt);
%         otherwise
%             error('GSP_DESIGN_WARPED_TRANSLATES: Unknown base filter')
%     end
% 
%     % Warp the uniform translates by log to generate the "wavelet" filters
%     wavelet_filters=cell(Nf-1,1);
%     for i=1:Nf-1
%         wavelet_filters{i} = @(x) uniform_filters{i}(wflog(x)); 
%     end
% %     g=cell(Nf,1);
% %     g(2:Nf)=wavelet_filters;
% %     % Generate a "scaling" filter that results in a tight frame
% %     g{1}= gsp_tighten_filter(lmax, wavelet_filters );
%     g{2:end} = wavelet_filters;
% else
        


    % Generate uniform translates covering [0,param.upper_bound_translates]
    if ~iscell(param.filter)
        switch param.filter
            case 'half_cosine'
                uniform_filters = gsp_design_half_cosine(wf(lmax),Nf);
            case 'itersine'
                paramt.overlap = param.overlap;
                uniform_filters = gsp_design_itersine(wf(lmax),Nf,paramt);
            otherwise
                error('GSP_DESIGN_WARPED_TRANSLATES: Unknown base filter')
        end
    else
        uniform_filters = cell(Nf,1);
        for ii = 1:Nf
            uniform_filters{ii} = @(x) param.filter{ii}(x*lmax/wf(lmax));
        end
    end
    
    g=cell(Nf,1);

    % Generate custom warped filters (may be adapated to the specific spectrum via the custom warp function)
    for i=1:Nf
        g{i} = @(x) uniform_filters{i}(wf(x)); 
    end
    
% end


end
