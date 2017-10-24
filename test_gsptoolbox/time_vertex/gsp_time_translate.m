function [ ft ] = gsp_time_translate(G, f,tau,param)
%GSP_TIME_TRANSLATE Generalized time-vertex translation of the signal f to the node i and time t
%   Usage: ft = gsp_time_translate(G, f,tau,param);
%
%   Input parameters
%       G   : Time-Vertex Graph structure
%       f   : Time-Vertex signal
%       tau : Time location
%   Output parameters
%       ft  : Translated signal
%
%   This function translate the time-vertex signal *f* at time tau.
%
%   Additional parameters
%   ---------------------
%   * *param.boundary*  : Time boundary condition for the translation: 
%                                                'periodic' (fft),  (default)
%                                                'reflecting' (dct),
%                                                'absorbing' (zero-padding).
%
%
% Author :  Francesco Grassi

if nargin<5
    param=struct;
end

if ~isfield(param,'boundary'), param.boundary='periodic'; end

[N,T] = size(f);

error('Change the boundary!')

switch param.boundary
    case 'periodic'
        delta = gsp_delta(G.jtv.T,tau).';
        fhat = fft(f,[],2);
        operator = repmat(fft(delta,[],2),N,1);
        ft = sqrt(G.N)*ifft(fhat .* operator,[],2);
        
    case 'symmetric'
        delta = gsp_delta(G.jtv.T,tau).';
        fhat = dct(f.').';
        operator = repmat(dct(delta.').',N,1);
        ft = sqrt(G.N)*idct( (fhat .* operator).' ).';
    
    case 'absorbing'
        
        ind = round(-tau*G.jtv.fs);
        ft = circshift(f,[0 ind]);
        
        if tau >0
           
            ft = [ft(:,1:end+ind) zeros(N,-ind)];
        
        else
            
            ft = [zeros(N,ind-1) ft(:,ind+1:end)];
            
        end
        
    otherwise
        error('Unknown boundary condition');
        
end

end