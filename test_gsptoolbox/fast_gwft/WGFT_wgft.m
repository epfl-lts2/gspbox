function [ wgft_coefficients ] = WGFT_wgft(g,V,signal,varargin)

% WGFT coefficients are output in an NxN matrix. The ith column corresponds
% to translation to vertex i. The kth row corresponds to frequency
% \lambda_{k-1}

N=size(V,1);
wgft_coefficients=zeros(N,N);

if ~isempty(varargin)
    param=varargin{1};
    for i=1:N
        for k=1:N
            g_ik=WGFT_atom(g,V,i,k,param);
            wgft_coefficients(k,i)=g_ik'*signal;
        end
    end
else
    for i=1:N
        for k=1:N
            g_ik=WGFT_atom(g,V,i,k);
            wgft_coefficients(k,i)=g_ik'*signal;
        end
    end    
end


end

