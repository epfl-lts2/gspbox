function [ errors ] = test_graphs(  )
%TEST_GRAPHS This function test all the graphs

errors = 0;

% Test creation of graph
errors = errors + test_swiss_roll();
errors = errors + test_david();
errors = errors + test_ring();
errors = errors + test_path();
errors = errors + test_airfoil();
errors = errors + test_comet();
errors = errors + test_erdos_renyi();
errors = errors + test_minnesota();
errors = errors + test_low_stretch_tree();
errors = errors + test_random_regular();
errors = errors + test_sensor();
errors = errors + test_community();

errors = errors + test_2dgrid();
errors = errors + test_torus();
errors = errors + test_gsp_graph();



end

function bool = test_symetry(W)
    bool = sum(sum(abs(W-W')))>1e-8;
end




function error = test_gsp_graph()

error = 0;
try
   figure(100);
   G = gsp_logo();
   gsp_plot_graph(G);
   close(100);
   fprintf('GRAPH: LOGO  1 ok\n');
catch
    error = 1;
    warning('GRAPH: Error in the LOGO 1 test')
end

if test_symetry(G.W)
    error = error + 1;
    warning('GRAPH: LOGO : not symetric')
else
    fprintf('GRAPH: LOGO symetric\n');
end

end

function error = test_community()

error = 0;
try
   figure(100);
   G = gsp_community(500);
   gsp_plot_graph(G);
   close(100);
   fprintf('GRAPH: community  1 ok\n');
catch
    error = 1;
    warning('GRAPH: Error in the community 1 test')
end

if test_symetry(G.W)
    error = error + 1;
    warning('GRAPH: community : not symetric')
else
    fprintf('GRAPH: community symetric\n');
end

end



function error = test_swiss_roll()

error = 0;
try
   figure(100);
   G = gsp_swiss_roll(200);
   gsp_plot_graph(G);
   close(100);
   fprintf('GRAPH: swiss roll 1 ok\n');
catch
    error = 1;
    warning('GRAPH: Error in the Swiss roll 2 test')
end

if test_symetry(G.W)
    error = error + 1;
    warning('GRAPH: Swiss roll : not symetric')
else
    fprintf('GRAPH: swiss roll symetric\n');
end

end



function error = test_david()

error = 0;
try
   figure(100);
   G = gsp_david_sensor_network(64);
   gsp_plot_graph(G);
   close(100);
   fprintf('GRAPH: David 1 ok\n');
catch
    error = error + 1;
    warning('GRAPH: Error in the 1 David test')
end

   
   G = gsp_david_sensor_network(64);
   if gsp_check_connectivity(G)
        fprintf('GRAPH: David 2 ok\n');
   else
        error = error + 1;
        warning('GRAPH: Error in the 2 David test')
   end

   G = gsp_david_sensor_network(500);
   if gsp_check_connectivity(G)
        fprintf('GRAPH: David 3 ok\n');
   else
        error = error + 1;
        warning('GRAPH: Error in the 3 David test')
   end

   
if test_symetry(G.W)
    error = error + 1;
    warning('GRAPH: david not symetric')
else
    fprintf('GRAPH: david symetric\n');
end
   
end

function error = test_ring()

N = 64;
error = 0;
try
   figure(100);
   G = gsp_ring(N);
   gsp_plot_graph(G);
   close(100);
   fprintf('GRAPH: ring 1 ok\n');
catch
    error = 1;
    warning('GRAPH: Error ring 1 test')
end

G = gsp_compute_fourier_basis(G);

x = rand(N,1);

fftx = 1/sqrt(N)*fft(x);
hatx = gsp_gft(G,x);
fftx = sort(fftx);
hatx = sort(hatx);
if norm(fftx-hatx)<1e-12
       fprintf('GRAPH: ring 2 ok\n');
       
else
    error = error +1;
    warning('GRAPH: Error ring 2 test')
    norm(gsp_gft(G,x)-fftx(inds))
    
end


if test_symetry(G.W)
    error = error + 1;
    warning('GRAPH: ring not symetric')
else
    fprintf('GRAPH: ring symetric\n');
end

end


function error = test_path()

N = 64;
error = 0;
try
   figure(100);
   G = gsp_path(N);
   gsp_plot_graph(G);
   close(100);
   fprintf('GRAPH: path 1 ok\n');
catch
    error = 1;
    warning('GRAPH: Error path 1 test')
end

G = gsp_compute_fourier_basis(G);

x = rand(N,1);

fftx = dct(x);
hatx = gsp_gft(G,x);
% fftx = sort(fftx);
% hatx = sort(hatx);
if norm(fftx-hatx)<1e-12
       fprintf('GRAPH: path 2 ok\n');
       
else
    error = error +1;
    warning('GRAPH: Error path 2 test')
    norm(gsp_gft(G,x)-fftx(inds))
    
end


if test_symetry(G.W)
    error = error + 1;
    warning('GRAPH: path not symetric')
else
    fprintf('GRAPH: path  symetric\n');
end

end

function error = test_airfoil()

error = 0;
try
   figure(100);
   G = gsp_airfoil();
   gsp_plot_graph(G);
   close(100);
   fprintf('GRAPH: aifoil ok\n');
catch
    error = 1;
    warning('GRAPH: Error in the airfoil test')
end


if test_symetry(G.W)
    error = error + 1;
    warning('GRAPH: airfoil not symetric')
else
    fprintf('GRAPH: airfoil symetric\n');
end

end


function error = test_comet()

error = 0;
try
   figure(100);
   G = gsp_comet(64,14);
   gsp_plot_graph(G);
   close(100);
   fprintf('GRAPH: comet ok\n');
catch
    error = 1;
    warning('GRAPH: Error in the comet test')
end


if test_symetry(G.W)
    error = error + 1;
    warning('GRAPH: comet not symetric')
else
    fprintf('GRAPH: comet symetric\n');
end

end

function error = test_erdos_renyi()

error = 0;
try
   G = gsp_erdos_renyi(100,0.05);

   fprintf('GRAPH: erdos renyi ok\n');
catch
    error = 1;
    warning('GRAPH: Error in erdos renyi test')
end


if test_symetry(G.W)
    error = error + 1;
    warning('GRAPH: erdos renyi not symetric')
else
    fprintf('GRAPH: erdos renyi symetric\n');
end

end

function error = test_minnesota()

error = 0;
try
   figure(100);
   G = gsp_minnesota(1);
   gsp_plot_graph(G);
   close(100);
   fprintf('GRAPH: minnesota ok\n');
catch
    error = 1;
    warning('GRAPH: Error in minnesota test')
end


if test_symetry(G.W)
    error = error + 1;
    warning('GRAPH: minnesota not symetric')
else
    fprintf('GRAPH: minnesota symetric\n');
end


end

function error = test_low_stretch_tree()

error = 0;
try
   figure(100);
   G = gsp_low_stretch_tree(8);
   gsp_plot_graph(G);
   close(100);
   fprintf('GRAPH: low strech tree ok\n');
catch
    error = 1;
    warning('GRAPH: Error in low strech tree test')
end


if test_symetry(G.W)
    error = error + 1;
    warning('GRAPH: low strecht tree not symetric')
else
    fprintf('GRAPH: low strech tree symetric\n');
end

end

function error = test_random_regular()

error = 0;
try
   G = gsp_random_regular(100,4);

   fprintf('GRAPH: random regular ok\n');
catch
    error = 1;
    warning('GRAPH: Error in random regular test')
end


if test_symetry(G.W)
    error = error + 1;
    warning('GRAPH: random regular not symetric')
else
    fprintf('GRAPH: random regular symetric\n');
end

end

function error = test_sensor()

error = 0;
try
   figure(100);
   G = gsp_sensor(50);
   param.show_edges = 1;
   gsp_plot_graph(G,param);
   close(100);
   fprintf('GRAPH: sensor ok\n');
catch
    error = 1;
    warning('GRAPH: Error in sensor test')
end


if test_symetry(G.W)
    error = error + 1;
    warning('GRAPH: sensor not symetric')
else
    fprintf('GRAPH: sensor symetric\n');
end

end



function error = test_2dgrid()

N = 4;
M = 8;
error = 0;
try
   figure(100);
   G = gsp_2dgrid();
   G = gsp_2dgrid(N);
   G = gsp_2dgrid(N,M);
   G = gsp_2dgrid(M,N);
   gsp_plot_graph(G);
   close(100);
   fprintf('GRAPH: 2d grid 1 ok\n');
catch
    error = 1;
    warning('GRAPH: Error 2d grid 1 test')
end
% 
% G = gsp_compute_fourier_basis(G);
% 
% x = rand(N,M);
% 
% fftx = dct2(x);
% hatx = gsp_gft(G,x(:));
% fftx = sort(fftx(:));
% hatx = sort(hatx);
% if norm(fftx(:)-hatx)<1e-12
%        fprintf('GRAPH: 2d grid 2 ok\n');
%        
% else
%     error = error +1;
%     warning('GRAPH: Error 2d grid 2 test')
%     norm(gsp_gft(G,x(:))-fftx(:))
%     
% end


if test_symetry(G.W)
    error = error + 1;
    warning('GRAPH: 2d grid not symetric')
else
    fprintf('GRAPH: 2d grid symetric\n');
end

end



function error = test_torus()

N = 4;
M = 3;
error = 0;
try
   figure(100);
   G = gsp_torus();
   G = gsp_torus(N);
   G = gsp_torus(N,M);
   G = gsp_torus(M,N);
   gsp_plot_graph(G);
   close(100);
   fprintf('GRAPH: torus 1 ok\n');
catch
    error = 1;
    warning('GRAPH: Error torus 1 test')
end



if test_symetry(G.W)
    error = error + 1;
    warning('GRAPH: torus not symetric')
else
    fprintf('GRAPH: torus symetric\n');
end

% 
% G = gsp_compute_fourier_basis(G);
% 
% x = rand(N,M);
% 
% fftx = dct2(x);
% hatx = gsp_gft(G,x(:));
% fftx = sort(real(fftx(:)))./sqrt(M*N);
% hatx = sort(hatx);
% if norm(fftx(:)-hatx)<1e-12
%        fprintf('GRAPH: 2d ring 2 ok\n');
%        
% else
%     error = error +1;
%     warning('GRAPH: Error 2d ring 2 test')
%     norm(gsp_gft(G,x(:))-fftx(:))
%     
% end

end
