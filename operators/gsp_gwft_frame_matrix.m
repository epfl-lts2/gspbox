function [ F ] = gsp_gwft_frame_matrix(G,g,param )
%GSP_GWFT_FRAME_MATRIX Create the matrix of the GWFT frame
%   Usage:  F = gsp_fast_gwft(G, g param );
%           F = gsp_fast_gwft(G,g);
%
%   Input parameters:
%         G     : Graph
%         g     : window
%         param : Structure of optional parameter
%   Output parameters:
%         F     : Frame
%
%   This function compute the graph windowed Fourier transform frame. In
%   might require a lot of memory. we create a matrix of size N^2*N. 
%
%   *param* a Matlab structure containing the following fields:
% 
%     * *param.verbose* : 0 no log, 1 print main steps, 2 print all steps.
%     By default, it is 1.



%   AUTHOR : Nathanael Perraudin
%   TESTING: test_fast_gwft

% Optional input arguments
if nargin<4, param=struct; end
if ~isfield(param, 'verbose'), param.verbose=1 ; end

if param.verbose && G.N>256
    warning('Create a big matrix, you can use other methods.');
end

% Compute the translate of g
ghat=G.U'*g;
Ftrans=sqrt(G.N)*G.U*(repmat(ghat,1,G.N).*G.U');

F=gsp_repmatline(Ftrans,1,G.N) .* ...
    repmat(repmat(1./G.U(:,1),1,G.N).*G.U,1,G.N);


end

