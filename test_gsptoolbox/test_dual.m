function [ errors ] = test_dual( )
%TEST_DUAL 

errors = 0;

errors = errors+ test_duality();

errors = errors+ test_approx_filter();
errors = errors+ test_exact();
errors = errors+ test_evaluate();

end

function [errors] = test_duality()
    errors = 0;
    
    N = 100;
    G = gsp_sensor(N);
    G = gsp_estimate_lmax(G);
    g = gsp_design_abspline(G,8);
    gd = gsp_design_can_dual(g);
    
    
if gsp_test_duality(G, g,gd )
   fprintf('DUAL: duality ok\n');
else
    errors = errors + 1;
    warning('DUAL: error in test duality')
end

end

function [errors] = test_approx_filter()
    errors = 0;
    
    N = 100;
    order = 15;
    G = gsp_sensor(N);
    G = gsp_estimate_lmax(G);
    g = gsp_design_abspline(G,8);
    ga = gsp_approx_filter(G,g,order);
    paramplot.show_sum = 0;
    figure(1)
    gsp_plot_filter(G,g,paramplot);
    title('Original filters')
    figure(2)
    gsp_plot_filter(G,ga,paramplot);
    title('Approximate filters');

    x = rand(N,1);
    param.order = order;
    c1 = gsp_filter_analysis(G,g,x,param);
    c2 = gsp_filter_analysis(G,ga,x,param);
    close all;
    
if norm(c1-c2)/norm(c1) < eps(1000)
   fprintf('DUAL: approx filter ok\n');
else
    errors = errors + 1;
    norm(c1-c2)/norm(c1)
    warning('DUAL: error in approx filter test')
end

end

function [errors] = test_exact()
    errors = 0;
    
    N = 100;
    G = gsp_sensor(N);
    G = gsp_compute_fourier_basis(G);
    g = gsp_design_abspline(G,8);
    gd = gsp_design_can_dual(g);
    paramplot.show_sum = 0;
    figure(1)
    gsp_plot_filter(G,g,paramplot);
    title('Original filters')
    figure(2)
    gsp_plot_filter(G,gd,paramplot);
    title('Canonical dual filters');

    x = rand(N,1);
    param.method = 'exact';
    coeff = gsp_filter_analysis(G,g,x,param);
    xs = gsp_filter_synthesis(G,gd,coeff,param);
    
    close all;
    
if norm(xs-x)/norm(x) < eps(1000)
   fprintf('DUAL: exact ok\n');
else
    errors = errors + 1;
    norm(xs-x)/norm(x)
    warning('DUAL: error in exact test')
end

end


function [errors] = test_evaluate()
    errors = 0;
    Nx = 1000;
    N = 100;
    G = gsp_sensor(N);
    G = gsp_estimate_lmax(G);
    g = gsp_design_abspline(G,8);
    
    x = rand(Nx,1);
    
    M = length(g);
    %t1 = tic;
    % Compute coefficient of g
    gcoeff = gsp_filter_evaluate(g,x);
    %toc(t1);
    % Compute coefficient of h
    %t2 = tic;
    x2 = zeros(Nx,M);
    for ii = 1:Nx
        x2(ii,:) =  pinv(gcoeff(ii,:)'); 
    end
    %toc(t2)
    x1 = gsp_evaluate_can_dual( g,x );
    
    if norm(x2-x1,'fro')/norm(x1,'fro') < eps(1000)
        fprintf('DUAL: evaluate ok\n');
    else
        errors = errors + 1;
        norm(x2-x1,'fro')/norm(x1,'fro')
        warning('DUAL: error in evaluate test')
    end

end
