function [ errors ] = test_sparsify(  )
    gsp_reset_seed(0);

errors = 0;
errors = errors + test1();

errors = errors + test2();

end

function [ errors ] = test1( )

errors = 0;
try
    epsilon = 0.4;
    param.distribute = 1;
    param.Nc = 20;
    G = gsp_sensor(256,param);
    G2 = gsp_graph_sparsify(G,epsilon);
    figure(100);
    gsp_plot_graph(G);
    title('Original graph')
    figure(101);
    gsp_plot_graph(G2);
    title('Sparsified graph')
    close(100);
    close(101);
   fprintf('SPARSIFY: test 1 ok\n');
catch
    errors = errors +1;
    warning('SPARSIFY: test 1 error')
end

end


function [ errors ] = test2( )

errors = 0;
    N = 100;
    epsilon = 0.6;
    G = gsp_sensor(N);
    gsp_reset_seed(0);
    G2 = gsp_graph_sparsify(G,epsilon);
    gsp_reset_seed(0);
    L = gsp_graph_sparsify_old(G.L,epsilon);

if sum(sum(abs(G2.L-L)))<1e-10
    
   fprintf('SPARSIFY: test  ok\n');
else
    errors = errors+1;
    badness = sum(sum(abs(G2.L-L)))
    warning('SPARSIFY: test 2 error')
end

end
