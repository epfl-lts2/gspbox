function [ errors ] = gsp_test_tig(  )

errors = 0;


Nf = 5;
N = 10;
G = gsp_sensor(N);
G = gsp_estimate_lmax(G);
g = gsp_design_itersine(G,Nf);


ntig1 = gsp_norm_tig( G,g, 1 );
ntig2 = gsp_norm_tig( G,g, 0,50000 );

if norm(ntig1(:)-ntig2(:))/norm(ntig1)<1e-2
   fprintf('TIG test 1 ok\n')
else
   norm(ntig1(:)-ntig2(:))/norm(ntig1)
   errors = errors +1;
   fprintf('Error in TIG test 1 !!!!!!!!!!\n')
end

ntig3 = gsp_norm_tig( G,g, 1,1 );
ntig4 = gsp_norm_tig( G,g, 1,2 );
ntig5 = gsp_norm_tig( G,g, 1,3 );
ntig6 = gsp_norm_tig( G,g, 1,-1 );



if norm(ntig1(:)-ntig3(:))/norm(ntig1)<1e-10
   fprintf('TIG test 2 ok\n')
else
   norm(ntig1(:)-ntig3(:))/norm(ntig1)
   errors = errors +1;
   fprintf('Error in TIG test 2 !!!!!!!!!!\n')
end

if norm(ntig1(:)-ntig4(:))/norm(ntig1)<1e-10
   fprintf('TIG test 3 ok\n')
else
   norm(ntig1(:)-ntig4(:))/norm(ntig1)
   errors = errors +1;
   fprintf('Error in TIG test 3 !!!!!!!!!!\n')
end

if norm(ntig1(:)-ntig5(:))/norm(ntig1)<1e-10
   fprintf('TIG test 4 ok\n')
else
   norm(ntig1(:)-ntig5(:))/norm(ntig1)
   errors = errors +1;
   fprintf('Error in TIG test 4 !!!!!!!!!!\n')
end


if norm(ntig1(:)-ntig6(:))/norm(ntig1)<1e-10
   fprintf('TIG test 5 ok\n')
else
   norm(ntig1(:)-ntig6(:))/norm(ntig1)
   errors = errors +1;
   fprintf('Error in TIG test 5 !!!!!!!!!!\n')
end

end

