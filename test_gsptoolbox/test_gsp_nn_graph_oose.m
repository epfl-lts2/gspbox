function [ errors ] = test_gsp_nn_graph_oose(  )

errors = 0;


Nx =10;
N = 30;
a = 0.1;
%% Random sampling
xg = rand(N,2);

%%
vx = linspace(-0.2,1.2,Nx);
vy = linspace(-0.2,1.2,Nx);
[X, Y] = meshgrid(vx, vy );
xo = [X(:),Y(:)];

%% Create the graph
param.sigma = 0.01;
param.epsilon = sqrt(-log(a)*param.sigma);
param.type = 'radius';
[ G ] = gsp_nn_graph_oose( xg, xo , param );
 t = param.sigma; 
G2  = gsp_create_cont_expW( xg,xo,t,a);

d1 = G.Gm.W(:);
d2 = G2.Gm.W(:);

error('This test makes no sense for now')

if norm(d1-d2,'fro')/norm(d1)<1e-8
    fprintf('TEST GSP_NN_GRAPH_OOSE OK\n');
    
else
    fprintf('ERROR IN TEST GSP_NN_GRAPH_OOSE !!!!!!!\n')
    norm(d1-d2,'fro')/norm(d1)  
    errors = errors +1;
end

end

