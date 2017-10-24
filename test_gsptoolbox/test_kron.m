function [ errors ] = test_kron(  )

errors = 0;

errors = errors + test_self_loop();
errors = errors + test_connected();

errors = errors + test_result();

end

function [ errors ] = test_self_loop( )

tol = 1e-10;

errors = 0;
N = 100;
ind = 1:2:N;
G = gsp_sensor(N);
L = gsp_kron_reduction(G.L,ind);
W = diag(diag(L)) - L;
D = sum(W);

if norm(full(D')-full(diag(L)))<tol
    fprintf('KRON: test self loop ok\n');
else
    warning('KRON: error in test self loop');
    errors = errors+1;
end
end

function [ errors ] = test_connected( )


errors = 0;
N = 100;
ind = 1:2:N;
G = gsp_sensor(N);
G = gsp_kron_reduction(G,ind);

if gsp_check_connectivity(G);
    fprintf('KRON: test connectivity ok\n');
else
    warning('KRON: error test connectivity');
    errors = errors+1;
end

end


function [ errors ] = test_result( )


errors = 0;
N = 100;
ind = 1:2:N;
G = gsp_sensor(N);
Gn = gsp_kron_reduction(G,ind);
L = gsp_kron_reduce_old(G.L,ind);


if sum(sum(abs(Gn.L-L)))<1e-12;
    fprintf('KRON: test result old ok\n');
else
    warning('KRON: test result old');
    badness = sum(sum(abs(Gn.L-L)))
    errors = errors+1;
end

end
