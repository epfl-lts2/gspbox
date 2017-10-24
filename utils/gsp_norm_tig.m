function [ ntig ] = gsp_norm_tig( G,g, exact,M,param )
%GSP_NORM_TIG Compute the norm of the tig of a frame
%   Usage:  gsp_norm_tig( G,g );
%           gsp_norm_tig( G,g, exact );
%
%   Input parameters:
%       G       : Graph
%       g       : filterbank
%       exact   : Exact method (default 1)
%       M       : Order for the approximation (default 50)
%   Output parameters:
%       ntig    : norm of tig
%
%   This function compute the norm of all atoms of the filterbank g.
% 
%   If *exact* is set to one, you can compute the norm chunk by chunk by
%   setting *M* to the size of the desired chunk. If *M* is -1, then
%   parpool is used.
%

% Author: Nathanael Perraudin
% Date  : 22 December 2014
% Tesing: gsp_test_tig

if nargin<3
    exact = 1;
end

if nargin<4
    if exact
        M = G.N; 
    else
        M = 50;
    end
end

if nargin<5
    param = struct;
end


if iscell(G)
    NG = length(G);
    ntig = cell(NG,1);
    for ii =1:NG
        ntig{ii} = gsp_norm_tig( G{ii},g{ii}, exact,M ,param);
    end
    return
end


Nf = length(g);
N = G.N;

if exact
    if M == -1
        if isempty(gcp('nocreate'))
            parpool
        end
%         parfor n = 1:N;
%             xin = zeros(N,1);
%             xin(n)= 1;
%             tig = gsp_vec2mat(gsp_filter_analysis(G,g,xin,param),Nf);
%             ntigt = sum(abs(tig).^2);
%             ntig(:,n) = reshape(ntigt,Nf,1);
%         end
        M = abs(M);
        N = ceil(G.N / M); % Chunks of M to save runtime memory.
        ntig = zeros(Nf,G.N);
        ntig2 = zeros(Nf*M,(N-1));
        parfor n = 1:(N-1);
            xin = zeros(G.N,M);
            for ii = 1:M
                xin(ii+(n-1)*M,ii)= 1;
            end
            tig = gsp_vec2mat(gsp_filter_analysis(G,g,xin),Nf);
            ntigt = sum(abs(tig).^2);
            ntig2(:,n) = reshape(ntigt,Nf*M,1);        
        end
        for n = 1:(N-1)
            ntig(:,(1:M)+(n-1)*M) = reshape(ntig2(:,n),Nf,M);            
        end
        M = G.N - M * (N-1);
        xin = zeros(G.N,M);
        for ii = 1:M
            xin(ii+G.N-M,ii) = 1;
        end
        tig = gsp_vec2mat(gsp_filter_analysis(G,g,xin),Nf);
        ntigt = sum(abs(tig).^2);
        ntig(:,(end-M+1) : end) = reshape(ntigt,Nf,M);
    else
        N = ceil(G.N / M); % Chunks of 1000 to save runtime memory.
        ntig = zeros(Nf,G.N);
        for n = 1:(N-1);
            xin = zeros(G.N,M);
            for ii = 1:M
                xin(ii+(n-1)*M,ii)= 1;
            end
            tig = gsp_vec2mat(gsp_filter_analysis(G,g,xin,param),Nf);
            ntigt = sum(abs(tig).^2);
            ntig(:,(1:M)+(n-1)*M) = reshape(ntigt,Nf,M);
        end
        M = G.N - M * (N-1);
        xin = zeros(G.N,M);
        for ii = 1:M
            xin(ii+G.N-M,ii) = 1;
        end
        tig = gsp_vec2mat(gsp_filter_analysis(G,g,xin,param),Nf);
        ntigt = sum(abs(tig).^2);
        ntig(:,(end-M+1) : end) = reshape(ntigt,Nf,M);
    end
    
    ntig = ntig';
else
    WN = sign(randn(N,M));
    filter_WN = gsp_vec2mat(gsp_filter_analysis(G,g,WN,param),Nf);
    ntig = sum(abs(filter_WN).^2,3)/M;
end

ntig = sqrt(ntig);
end

