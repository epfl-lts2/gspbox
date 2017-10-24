function [ C ] = gsp_gwft(G,g,f, param )
%GSP_GWFT Generalized windowed Fourier transform
%   Usage:  C = gsp_gwft(G,g,f, param );
%           C = gsp_gwft(G,g,f);
%
%   Input parameters:
%         G     : Graph
%         g     : Window (graph signal or kernel)
%         f     : Graph signal (column vector)
%         param : Structure of optional parameter
%   Output parameters:
%         C     : Coefficient.
%
%   This function compute the graph windowed Fourier transform of a signal
%   *f* with the window *g*. The function returns a matrix of size N^2*N. 
%
%   *param* a Matlab structure containing the following fields:
% 
%     * *param.verbose* : 0 no log, 1 print main steps, 2 print all steps.
%       By default, it is 1.
%     * *param.lowmemory* : use less memory. By default, it is 1.
%
%   Reference: shuman2013windowed

% Author : Nathanael Perraudin
% Testing: test_gwft


% Optional input arguments
if nargin<4, param=struct; end
if ~isfield(param, 'verbose'), param.verbose=1 ; end
if ~isfield(param, 'lowmemory'), param.lowmemory=1 ; end

if ~isfield(G,'U')
   error(['GSP_GWFT: You need first to compute the Fourier basis. ',...
       'You can do it with the function gsp_compute_fourier_basis']);
end

if sum(G.U(:,1)<eps(2))
    error(['GSP_GWFT: The current implementation of this function is ',...
       'not working for disconnected graphs'])
end

Nf = size(f,2);

if iscell(g)
    g = gsp_igft(G,g{1}(G.e));
end


if isa(g, 'function_handle')
    g = gsp_igft(G,g(G.e));
end

if ~param.lowmemory
    % Compute the Frame into a big matrix
    Frame=gsp_gwft_frame_matrix(G,g,param);

    C=Frame'*f;

    C=reshape(C,G.N,G.N,Nf);
else
    % Compute the translate of g
    ghat=G.U'*g;
    Ftrans=sqrt(G.N)*G.U*(repmat(ghat,1,G.N).*G.U');
    C=zeros(G.N);
    for jj = 1:Nf
        for ii=1:G.N
             C(:,ii,jj)=(repmat(1./G.U(:,1),1,G.N) .* ...
                      G.U.*repmat(Ftrans(:,ii),1,G.N))'*f(:,jj);
        end
    end
    
end

