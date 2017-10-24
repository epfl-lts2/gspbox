% GSPBOX - Operators
%
%  Localisation
%    gsp_localize       -  Localize a kernel
%    gsp_modulate       -  Generalized modulation operator
%
%  Differential
%    gsp_grad_mat       -  Compute the gradient sparse matrix
%    gsp_grad           -  Compute the gradient of a signal
%    gsp_div            -  Compute the divergence of a signal
%
%  Transforms
%    gsp_gft            -  Graph Fourier transform
%    gsp_igft           -  Inverse graph Fourier transform
%    gsp_gwft           -  Windowed graph Fourier transform
%    gsp_ngwft          -  Normalized windowed graph Fourier transform
%
%  Time-Vertex Transforms
%    gsp_jft            -  Joint Time-Vertex Fourier transform
%    gsp_ijft           -  Inverse Joint Time-Vertex Fourier transform
%    gsp_tft            -  Time-Vertex Time-Fourier transform
%    gsp_itft           -  Time-Vertex Inverse Time-Fourier transform
%
%  Pyramid - Reduction
%    gsp_kron_reduce    - Kron reduction
%    gsp_graph_multiresolution - Compute a multiresolution of graphs
%    gsp_pyramid_analysis - Analysis operator for graph pyramid
%    gsp_pyramid_analysis_single - Compute a single level of the graph pyramid transform coefficients
%    gsp_pyramid_synthesis - Sythesis operator for graph pyramid
%    gsp_pyramid_synthesis_single -Synthesize a single level of the graph pyramid transform 
%    gsp_pyramid_cell2coeff - Keep only the necessary coefficients
%    gsp_interpolate    - Interpolate a signal
%
%  For help, bug reports, suggestions etc. please send email to
%  gspbox 'dash' support 'at' groupes 'dot' epfl 'dot' ch
%
% see also: prox

