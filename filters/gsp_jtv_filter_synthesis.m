function X = gsp_jtv_filter_synthesis(G, g,filtertype, c, param)
%GSP_JTV_FILTER_SYNTHESIS Synthesis operator of a time-vertex filterbank
%   Usage:  X = gsp_jtv_filter_synthesis(G, g, c, param)
%
%   Input parameters:
%         G          : Time-Vertex graph structure
%         g          : Cell array of time-vertex filters
%         filtertype : Filter domain (ts,js,ts-array,js-array)
%         c          : Coefficients matrix
%         param      : Struct of parameters
%   Output parameters:
%         X          : Time-Vertex signal
%
%   Synthesis operator of a time-vertex filterbank
%
%   Additional parameters
%   ---------------------
%   * *param.verbose*     : verbosity level. 0 no log - 1 display warnings.
%     (default 1)
%   * *param.method*      : (default 'exact')
%   * *param.vectorize*   : 1 if coefficients are in vectorized form (default 0)
%

% Author :  Francesco Grassi
% Date : July 2016


% Read input parameters
if nargin < 5
    param = struct;
end


if ~gsp_check_filtertype(filtertype)
    error('Invalid filtertype');
end


if ~isfield(param,'method')
    if isfield(G,'U')
        param.method = 'exact';
    else
        param.method = 'cheby';
    end
end

if ~isfield(param,'order');     param.order = 40; end
if ~isfield(param,'verbose');   param.verbose = 1; end
if ~isfield(param,'vectorize'), param.vectorize = 0;end


if ~gsp_check_jtv(G)
    error('GSP_JTV_FILTER_SYNTHESIS need the time dimension. Use GSP_JTV_GRAPH.');
end

N   = G.N;
T   = G.jtv.T;
fs  = G.jtv.fs;
lag = size(c,2);

if isnumeric(g)
    Nf = size(g,3);
else
    switch filtertype
        case {'js-array','ts-array'}
            Nf = size(g,1);
        case {'ts', 'js'}
            Nf = numel(g);
    end
end

Ns  = size(c,4);

switch filtertype
    case  {'ts','ts-array'}
        t = 0:1/fs:(lag-1)/fs;
    case  {'js','js-array'}
        t = gsp_cfa(lag,fs);
    otherwise
        error('Unknown filtertype');
end

if param.vectorize
    c=gsp_jtv_vec2mat(G,c);
end
switch param.method
    
    case 'exact'
        
        
        
        param.domain = 'joint-spectral'; %output filter in joint spectral domain
        
        fie = gsp_jtv_filter_evaluate(g,filtertype,G.e,t,param);
        
        X = zeros(N,T,Ns);
        
        for n=1:Nf
            
            tmp = gsp_ijft(G,gsp_jft(G,c(:,:,n,:)).*repmat(fie(:,:,n),1,1,1,Ns));
            
            %if imaginary part is small remove it
            if sum(abs(imag(tmp(:))))<(1e-10 * norm(tmp(:),'fro'));
                tmp = real(tmp);
            end
            
            tmp = tmp(:,lag-T+1:end,:,:);
            
            X = X + reshape(tmp,G.N,T,Ns);
        end
        
    case 'cheby'
        
        X = zeros(N,T,Ns);
        cdot = gsp_tft(G,c);
        
        switch filtertype
            case {'ts','ts-array'}
                cheby_coeffs = gsp_jtv_cheby_coeff(G,g,filtertype,param.order, param.order +1);
                
                for n=1:Nf
                    X = X + gsp_jtv_cheby_op(G, cheby_coeffs(:,:,n), squeeze(cdot(:,:,n,:)));
                end
                
            case {'js','js-array'}
                
                for ii = 1:T
                    X(:,ii,:) = gsp_filter_synthesis(G,@(x) g(x,t(ii)),squeeze(cdot(:,ii,:,:)),param);
                end
                
                
        end
        
        X = gsp_itft(G,X);
        
end

if param.vectorize
    X=gsp_jtv_mat2vec(X);
end

end
