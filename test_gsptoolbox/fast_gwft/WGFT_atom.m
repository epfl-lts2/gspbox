function [ g_ik ] = WGFT_atom(g,V,i,k,varargin)

if ~isempty(varargin)
    param=varargin{1};
else
    param=0;
end

% If not computing eigenvectors, pass L in place of V - should change
% internal label here

N=size(V,1);

% if isfield(param, 'mod_first')
%     mod_first=param.mod_first;
% else
%     mod_first=0;
% end
% 
% if isfield(param, 'mod_abs')
%     mod_abs=param.mod_abs;
% else
%     mod_abs=0;
% end

if isfield(param, 'wgft_type')
    wgft_type=param.wgft_type;
else
    wgft_type='MkTi';
end


% if mod_first
%     if mod_abs
%         modulated_window_vertex=WGFT_gmod(g,k,V);
%         modulated_window_gsd=gft(modulated_window_vertex,V);
%         abs_mod_window=igft(abs(modulated_window_gsd),V);
%         g_ik=WGFT_gtrans(abs_mod_window,i,V);
%     else     
%         g_ik=WGFT_gtrans(WGFT_gmod(g,k,V),i,V);
%     end
% else
%     g_ik=WGFT_gmod(WGFT_gtrans(g,i,V),k,V);
% end


switch wgft_type
    
    case 'MkTi'
        g_ik=WGFT_gmod(WGFT_gtrans(g,i,V),k,V);
    
    case 'TiMk'
        g_ik=WGFT_gtrans(WGFT_gmod(g,k,V),i,V);
        
    case 'TiabsMk'
        modulated_window_vertex=WGFT_gmod(g,k,V);
        modulated_window_gsd=gft(modulated_window_vertex,V);
        abs_mod_window=igft(abs(modulated_window_gsd),V);
        g_ik=WGFT_gtrans(abs_mod_window,i,V);
        
    case 'Ticustom'
        % Not using g or V
        if ~isfield(param, 'custom_approx_mod_kernels')
            error('Custom approx_mod_kernels not specified');
        end
%         if ~isfield(param, 'tau')
%             error('tau parameter not provided');
%         else
%             tau=param.tau;
%         end
        approx_modulated_kernel=param.custom_approx_mod_kernels{k};
        inds=[1:N];
        approx_modulated_kernel_gsd=approx_modulated_kernel(inds');
        approx_modulated_kernel_vertex=igft(approx_modulated_kernel_gsd,V);
        g_ik=WGFT_gtrans(approx_modulated_kernel_vertex,i,V);
    
    case 'ApproxTransContinuousKernel'
        L=V;
        if ~isfield(param, 'continuous_kernels')
            error('Continuous approx_mod_kernels not specified');
        end
        if ~isfield(param, 'lmax')
            lmax=sgwt_rough_lmax(L);
        else
            lmax=param.lmax;
        end
        if ~isfield(param, 'order')
            order=25;
        else
            order=param.order;
        end
        g_ik= WGFT_approx_trans(param.continuous_kernels{k},i,L,lmax,order);
    otherwise
        warning('wgft_type not recognized. Using MkTi instead');
        g_ik=WGFT_gmod(WGFT_gtrans(g,i,V),k,V);


end

