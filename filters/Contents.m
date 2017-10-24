% GSPBOX - Filters
%
%  Design - N filter
%    gsp_design_mexican_hat       -  Design a mexican hat filterbank
%    gsp_design_abspline          -  Design a abspline filterbank
%    gsp_design_meyer             -  Design a Meyer filterbank (tight)
%    gsp_design_simple_tf         -  Design a simple tight frame filterbank
%    gsp_design_itersine          -  Design a itersine filterbank (tight)
%    gsp_design_half_cosine       -  Design a half cosine filterbank (tight)
%    gsp_design_warped_translates -  Design a filterbank with a warping function 
%
%  Design - 2 filter (LP - HP) tight filterbank
%    gsp_design_regular           -  Design 2 filter with the "regular" construction
%    gsp_design_held              -  Design 2 filter with the "Held" construction
%    gsp_design_simoncelli        -  Design 2 filter with the "Simoncelli" construction
%    gsp_design_papadakis         -  Design 2 filter with the "Papadakis" construction
%
%  Dual filterbank
%    gsp_design_can_dual          -  Design the canonical dual filterbank
%    gsp_evaluate_can_dual        -  Evaluate the canonical dual of a filterbank
%    gsp_test_duality             -  Test if 2 filterbanks are dual
%    gsp_test_duality_coefficient -  Test if 2 discrete filterbanks are dual
%
%  Low pass filters
%    gsp_design_heat              -  Design a heat kernel filter
%    gsp_design_expwin            -  Design a expwin filter
%    gsp_design_smooth_indicator  -  Design a smooth indicator function
%
%  Application
%    gsp_filter                   -  Shortcut to gsp_filter_analysis
%    gsp_filter_analysis          -  Analysis operator for filterbank
%    gsp_filter_synthesis         -  Synthesis operator for filterbank
%    gsp_filter_inverse           -  Inverse operator for filterbank
%
%  Joint Time-Vertex Filter Design
%    gsp_jtv_design_diffusion         -  Design a diffusion filterbank
%    gsp_jtv_design_wave              -  Design a wave filterbank
%    gsp_jtv_design_damped_wave       -  Design a damped wave filterbank
%    gsp_jtv_design_dgw               -  Design a generic dynamic graph wavelet
%
%  Joint Time-Vertex Filter Application
%    gsp_jtv_filter_analysis      -  Analysis operator for time-vertex filterbank
%    gsp_jtv_filter_synthesis     -  Synthesis operator for time-vertex filterbank
%    gsp_jtv_filter_evaluate      -  Evaluate a time-vertex filterbank
%    gsp_jtv_filter_array         -  Convert a ts/js time-vertex filter to a ts/js-array time-vertex filterbank
%    gsp_jtv_compute_frame        -  Return the matrix operator associated to a time-vertex filterbank
%    gsp_jtv_evaluate_can_dual    -  Evaluate the canonical dual of a time-vertex filterbank
%    gsp_jtv_design_can_dual      -  Design the canonical dual of a time-vertex filterbank
%    gsp_filter_inverse           -  Inverse operator for a time-vertex filterbank
%
%  Size Handling
%    gsp_mat2vec                  -  Matrix to vector representation for filterbanks 
%    gsp_vec2mat                  -  Vector to matrix representation for filterbanks 
%
%  Utils
%    gsp_approx_filter            -  Create an approximation of a filterbank with Chebyshev
%    gsp_wlog_scales              -  Compute log scale vector for wavelets
%    gsp_filter_evaluate          -  Evaluate a filterbank
%    gsp_filterbank_bounds        -  Bound for the filterbank
%    gsp_tighten_filter           -  Create a filter that tighten the filterbank
%    gsp_warp_filter              -  Warp a filter
%    gsp_multiply_filters         -  Multiply two filters
%    gsp_filterbank_matrix        -  Return the matrix operator associated to a filterbank
%
%  For help, bug reports, suggestions etc. please send email to
%  gspbox 'dash' support 'at' groupes 'dot' epfl 'dot' ch
%

% To be done
%   - add function in doc from meyer and simple_tf


