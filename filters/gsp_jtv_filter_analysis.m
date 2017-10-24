function [c] = gsp_jtv_filter_analysis(G, g, filtertype, X, param)
%GSP_JTV_FILTER_ANALYSIS Analysis operator of a time-vertex filterbank
%   Usage:  [c] = gsp_jtv_filter_analysis(G, g, X)
%           [c] = gsp_jtv_filter_analysis(G, g, X, param)
%
%   Input parameters:
%         G          : Time-Vertex graph structure
%         g          : Cell array of time-vertex filters
%         filtertype : Filter domain (ts,js,ts-array,js-array)
%         X          : Time-Vertex signal
%   Output parameters:
%         c          : Coefficients matrix
%
%   Analysis operator for a time-vertex filterbank
%
%   Additional parameters
%   ---------------------
%
%   * *param.verbose*   : verbosity level. 0 no log - 1 display warnings (default 1)
%   * *param.method*    : Exact or Cheby
%   * *param.order*     : Order of Cheby approximation
%   * *param.vectorize*   : 1 if coefficients are in vectorized form (default 0)
%


% Author :  Francesco Grassi
% Date : September 2016
% Testing: test_jtv_filter


% Read input parameters
if nargin < 5
    param = struct;
end


if ~isfield(param,'method')
    if isfield(G,'U')
        param.method = 'exact';
    else
        param.method = 'cheby';
    end
end

if ~isfield(param,'order'); param.order = 40; end
if ~isfield(param,'verbose'); param.verbose = 1; end
if ~isfield(param,'vectorize'), param.vectorize = 0;end

if ~gsp_check_jtv(G)
    error('GSP_JTV_FILTER_ANALYSIS need the time dimension. Use GSP_JTV_GRAPH.');
end

T = G.jtv.T;


switch filtertype
    case  {'ts','ts-array'}
        v = gsp_jtv_ta(G);
    case  {'js','js-array'}
        v = gsp_jtv_fa(G,0);
    otherwise
        error('Unknown filtertype');
end


if G.jtv.extension
    Ts = size(X,2);
    X = [ zeros(G.N,2*T-Ts-1) X ];
    zeropad = T-1;
    tau = 2*T-1;
else
    zeropad = 0;
    tau = T;
end

%Number of signals
Ns = size(X,3);

if param.vectorize
    X=gsp_jtv_vec2mat(G,X);
end

switch param.method
    case 'exact'
        
        param.domain = 'joint-spectral'; %output filter in joint spectral domain
        fie = gsp_jtv_filter_evaluate(g,filtertype,G.e,v,param);
        
        Nf  = size(fie,3);
        
        c = zeros(G.N,tau,Nf,Ns);
        
         if zeropad>0
            fie = gsp_jft(G,[gsp_ijft(G,fie) zeros(G.N,zeropad,Nf) ]);
         end 
        
        for n=1:Nf
            
                tmp = gsp_ijft(G,gsp_jft(G,X).*repmat(conj(fie(:,:,n)),1,1,Ns));
                
                %if imaginary part is small remove it
                if sum(abs(imag(tmp(:))))<(1e-10 * norm(tmp(:),'fro'));
                    tmp = real(tmp);
                end
                
                c(:,:,n,:) = reshape(tmp,G.N,tau,1,Ns);
            
        end
        
        
    case 'cheby'
        
        switch filtertype
            case {'ts','ts-array'}
                Xdot = gsp_tft(G,X);
                
                cheby_coeffs = gsp_jtv_cheby_coeff(G,g,filtertype,param.order, param.order +1);
                
                c = gsp_jtv_cheby_op(G, conj(cheby_coeffs), Xdot);
                c = gsp_itft(G,c);
                
                
                
            case {'js','js-array'}
                
                Xdot = gsp_tft(G,X);
                c    = zeros(G.N,T);
                
                for ii = 1:T
                    c(:,ii) = gsp_filter_analysis(G,@(x) g(x,v(ii)),Xdot(:,ii),param);
                end
                
                c = gsp_itft(G,c);
        end
        
end

if param.vectorize
    c = gsp_jtv_mat2vec(c);
end


end
