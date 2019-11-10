function [c] = gsp_filter_analysis(G, fi, s, param)
%GSP_FILTER_ANALYSIS Analysis operator of a gsp filterbank
%   Usage:  coeffs = gsp_filter_analysis(G, fi, signal);
%           coeffs = gsp_filter_analysis(G, fi, signal, param);
%
%   Input parameters:
%         G         : Graph structure.
%         fi        : Set of spectral graph filters.
%         s         : graph signal to analyze.
%         param     : Optional parameter
%   Output parameters:
%         c         : Transform coefficients
%
%   'gsp_filter_analysis(G,fi,signal)' computes the transform
%   coefficients of a signal $f$, where the atoms of the transform
%   dictionary are generalized translations of each graph spectral filter
%   to each vertex on the graph.
%
%   .. c = D' * f 
%
%   .. math:: c =  D^* f
%
%   where the columns of $D$ are $g_{i,m}=T_i g_m$, and $T_i$ is a
%   generalized translation operator applied to each filter 
%   $\hat{g}_m(\cdot)$.  
%
%   Each column of *c* is the response of the signal to one filter.
%
%   Example:::
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
%   * *param.method*  : Select the method to be used for the computation. 
%     * 'exact'     : Exact method using the graph Fourier matrix
%     * 'cheby'     : Chebyshev polynomial approximation
%     * 'lanczos'   : Lanczos approximation
%     Default: if the Fourier matrix is present: 'exact' otherwise 'cheby'
%   * *param.order* : Degree of the Chebyshev approximation
%     Default: 30
%   * *param.grid_order* : grid order used to compute quadrature 
%     Default: param.order+1
%   * *param.verbose* : Verbosity level (0 no log - 1 display warnings)
%     Default: 1   
%
%   See also: gsp_filter_synthesis gsp_filter_inverse 
% 
%   References: hammond2011wavelets
%

% Author: David I Shuman, Nathanael Perraudin
% Testing: test_filter
% Date: 19 March 2014

  
% Read input parameters
if nargin < 4
    param = struct;
end

if iscell(G)
    NG = numel(G);
    c = cell(NG,1);
    for ii = 1:NG
        warning('Check what happen here')
        if iscell(s)
            c{ii} = gsp_filter_analysis(G{ii}, fi{ii}, s{ii}, param);
        else
            c{ii} = gsp_filter_analysis(G{ii}, fi{ii}, s, param);
        end
    end
    return
end

if isnumeric(fi)
    Nf = size(fi,2);
else    
    Nf = numel(fi);
end

if isfield(param, 'exact')
    warning('param.exact is not used anymore. Please use param.method instead');
    if param.exact
        param.method = 'exact';
    else
        param.method = 'cheby';
    end
end

if ~isfield(param,'method')
    if gsp_check_fourier(G)
        param.method = 'exact';
    else
        param.method = 'cheby';
    end
end

if ~isfield(param,'order'); param.order = 30; end
if ~isfield(param,'grid_order'); param.grid_order = param.order+1; end
if ~isfield(param,'verbose'); param.verbose = 1; end

if isfield(param, 'cheb_order')
    param.order = param.cheb_order;
    warning('param.cheb_order is not used anymore. Please use param.order instead');
end



switch param.method
    case 'exact' 
        if ~gsp_check_fourier(G)
            if param.verbose
                warning(['GSP_FILTER_ANALYSIS: The Fourier matrix is not ',...
                    'available. The function will compute it for you. ',...
                    'However, if you apply many time this function, you ',...
                    'should precompute it using the function: ',...
                    'gsp_compute_fourier_basis']);
            end
            G = gsp_compute_fourier_basis(G);
        end
        Nv = size(s,2);
        c = zeros(G.N*Nf,Nv);
        if isnumeric(fi)
            fie = fi;
        else
            fie = gsp_filter_evaluate(fi,G.e);
        end
%         for ii=1:Nf
%             c((1:G.N)+G.N * (ii-1),:)= gsp_igft(G, ...
%                 repmat(conj(fie(:,ii)),1,Nv) ...
%                 .* gsp_gft(G, s));
%         end

%         c = gsp_mat2vec(gsp_igft(G, ...
%     repmat(conj(fie),[1,1,Nv]) ...
%     .* repmat(reshape(gsp_gft(G, s),G.N,1,[]),[1, Nf,1]) ));

shat = gsp_gft(G, s);
%chat = repmat(conj(fie),[1,1,Nv]) .* repmat(reshape(shat,G.N,1,[]),[1, Nf,1]);
chat = bsxfun(@times, conj(fie), permute(shat,[1 3 2]));

c = gsp_mat2vec(gsp_igft(G, chat));
        

    case 'cheby'

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
                param.order, param.grid_order, param);    
        c = gsp_cheby_op(G, cheb_coeffs, s);
%     case 'cheby_p'
% 
%         if ~isfield(G,'lmax');
%             G = gsp_estimate_lmax(G);
%             if param.verbose
%                 warning(['GSP_FILTER_ANALYSIS: The variable lmax is not ',...
%                     'available. The function will compute it for you. ',...
%                     'However, if you apply many time this function, you ',...
%                     'should precompute it using the function: ',...
%                     'gsp_estimate_lmax']);
%             end
%         end
% 
% 
%         cheb_coeffs = gsp_cheby_coeff(G, fi,...
%                 param.order, param.order +1);      
%         c = gsp_cheby_op_p(G, cheb_coeffs, s);    
    case 'lanczos'
        c = gsp_lanczos_op(G, fi, s, param);
   
    otherwise
        error('Unknown method: please select exact, cheby or lanczos');
end

end

