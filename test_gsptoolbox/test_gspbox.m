function [ errors ] = test_gspbox(  )
%RUN_TESTS this function run all test for the GSPBOX
%

clear all;
close all;
clc;

gsp_start();

errors=0;


errors = errors + test_resistance_distance();
errors = errors + test_graphs();
errors = errors + test_plotting();
errors = errors + test_operators();
errors = errors + test_filter();
errors = errors + test_gsp_prox();
errors = errors + test_gsp_dn();
errors = errors + test_sparsify();
errors = errors + test_kron();
errors = errors + test_pyramid();
errors = errors + test_gsp_solve_l1();
errors = errors + gsp_test_tig();
errors = errors + test_dual();
errors = errors + test_gsp_solve_l0();
errors = errors + test_gsp_solve_l1();
errors = errors + test_gsp_filter_manifold();
errors = errors + test_rmse();
errors = errors + test_gsp_hope_distanz();
errors = errors + test_graph_ml();
errors = errors + test_gsp_distanz();
% errors = errors + test_gsp_nn_graph_oose();
errors = errors + test_gsp_proj_filerbank();
errors = errors + test_complex();
errors = errors + test_wiener();
errors = errors + test_gsp_remove_mean();


errors = errors + test_jtv_filter();
errors = errors + test_gsp_wiener_jft_predition();
errors = errors + test_gsp_estimate_vertex_time_psd();






if errors
    fprintf('\n ****************************************** \n   They are %i test(s) ended with errors!\n ****************************************** \n',errors);
else
    fprintf('\n ****************************************** \n   No errors in testbenches!\n ****************************************** \n');
end

end

