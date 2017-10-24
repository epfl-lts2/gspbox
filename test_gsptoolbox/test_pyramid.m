function errors = test_pyramid()

gsp_reset_seed(0);
errors = 0;
errors = errors + test_landweber();

errors = errors + test_analysis();


errors = errors + test_result_graph();
errors = errors + test_result_graph2();
errors = errors + test_result_graph_multilevel();
errors = errors + test_result_graph_multilevel2();


errors = errors + test_interpolate();
errors = errors + test_synthesis();
errors = errors + test_reconstruction();


end

function errors = test_landweber()

errors = 0;
N = 100;
Nl = 3;
G = gsp_sensor(N);
lambda = 0.005;
g = @(x) 1./(1+x);
param.h_filters = {g,g,g};
param.lambda = lambda;
param.sparsify_epsilon = 0.5;
param.sparsity = 1;
Gs1 = gsp_graph_multiresolution(G,Nl,param);
Gs1 = gsp_estimate_lmax(Gs1);


% Test with 1d signal
f = rand(N,1);
f(1) = 10;
paramp.order = 100;
paramp.h_filters = param.h_filters;
paramp.landweber_its=100;
paramp.least_squares=1;

[ca1,pe1] = gsp_pyramid_analysis(Gs1,f,Nl, paramp);

% coeff = gsp_pyramid_cell2coeff(ca1,pe1);
s1 = gsp_pyramid_synthesis(Gs1,ca1{end},pe1,paramp);

if norm(s1-f) < 1e-10
    fprintf('PYRAMID: landweber ok\n');
else
    warning('PYRAMID: error in landweber');
    badness = norm(s1-f)
    errors = errors+1;
end

paramp.use_landweber = 0;
s1 = gsp_pyramid_synthesis(Gs1,ca1{end},pe1,paramp);

if norm(s1-f) < 1e-10
    fprintf('PYRAMID: least square full ok\n');
else
    warning('PYRAMID: error in least square full');
    badness = norm(s1-f)
    errors = errors+1;
end

end

function errors = test_synthesis()

errors = 0;

N = 10;
Nl = 1;
G = gsp_sensor(N);
lambda = 0.005;
g = @(x) 1./(1+x);
param.filters = g;
param.lambda = lambda;
param.sparsify_epsilon = 0.5;
param.sparsity = 1;
Gs1 = gsp_graph_multiresolution(G,Nl,param);
Gs1 = gsp_estimate_lmax(Gs1);
Gs2 = gsp_graph_multiresolution_old(G,Nl,param);
Gs2 = gsp_estimate_lmax(Gs2);

% Test with 1d signal
f = rand(N,1);
f(1) = 10;
paramp.order = 100;

[ca1,pe1] = gsp_pyramid_analysis(Gs1,f,Nl, paramp);
paramp.h_filters = g;
paramp.regularize_epsilon = lambda;
[ca2,pe2] = gsp_pyramid_analysis_old(f,Gs2,Nl,paramp);

% coeff = gsp_pyramid_cell2coeff(ca1,pe1);
s1 = gsp_pyramid_synthesis(Gs1,ca1{end},pe1,paramp);

s2 = gsp_pyramid_synthesis_old(ca2{end},pe2(1),Gs2,paramp);

if norm(s1-s2) < 1e-10
    fprintf('PYRAMID: test syntehis ok\n');
else
    warning('PYRAMID: error in test syntehis');
    badness = norm(s1-s2)
    errors = errors+1;
end

end

function errors = test_analysis()

errors = 0;

N = 10;
Nl = 1;
G = gsp_sensor(N);
lambda = 0.005;
g = @(x) 1./(1+x);
param.filters = g;
param.lambda = lambda;
Gs1 = gsp_graph_multiresolution(G,Nl,param);
Gs1 = gsp_estimate_lmax(Gs1);
Gs2 = gsp_graph_multiresolution_old(G,Nl,param);
Gs2 = gsp_estimate_lmax(Gs2);

% Test with 1d signal
f = rand(N,1);
f(1) = 10;
paramp.order = 100;
paramp.h_filters = g;
paramp.regularize_epsilon = lambda;
[ca1,pe1] = gsp_pyramid_analysis(Gs1,f, Nl,paramp);
[ca2,pe2] = gsp_pyramid_analysis_old(f,Gs2,Nl,paramp);



if sum(abs(ca1{2}-ca2{2})) < 1e-10
    fprintf('PYRAMID: test analysis 1d-1 ok\n');
else
    warning('PYRAMID: error in test  test analysis 1d-1');
    badness = sum(abs(ca1{2}-ca2{2})) 
    errors = errors+1;
end

if sum(abs(pe1{1}-pe2{1})) < 1e-10
    fprintf('PYRAMID: test analysis 1d-2 ok\n');
else
    warning('PYRAMID: error in test  test analysis 1d-2');
    badness = sum(abs(pe1{1}-pe2{1}))
    errors = errors+1;
end

% Test with 2d signal.
f2 = rand(N,3);
f2(1, :) = 10;
paramp.order = 100;

[ca3,pe3] = gsp_pyramid_analysis(Gs1,f2, Nl,paramp);
ca4 = cell(Nl+1, 1);
pa4 = cell(Nl+1, 1);
for ii = 1:3
    [tmpca4, tmppe4] = gsp_pyramid_analysis(Gs1,f2(:, ii), Nl,paramp);
    for jj =1:(Nl + 1);
        ca4{jj}(:, ii) = tmpca4{jj};
    end
    for jj =1:(Nl);
        pe4{jj}(:, ii) = tmppe4{jj};
    end

end

for ii = 1:Nl+1
    if ca3{ii} == ca4{ii};
        fprintf('PYRAMID: test analysis 2d-1 ok\n');
    else
        warning('PYRAMID: error in test  test analysis  2d-1');
        badness = sum(ca3{ii}-ca4{ii}) 
        errors = errors+1;
    end


end
for ii = 1:Nl
    if pe3{ii} == pe4{ii};
        fprintf('PYRAMID: test analysis 2d-2 ok\n');
    else
        warning('PYRAMID: error in test  test analysis 2d-2');
        badness = sum(abs(pe3{ii}-pe4{ii}))
        errors = errors+1;
    end
end
end

function errors = test_interpolate()

errors = 0;

N = 10;
Nl = 1;
G = gsp_sensor(N);
G = gsp_estimate_lmax(G);
lambda = 0.005;
param.lambda = lambda;
Gs = gsp_graph_multiresolution(G,Nl,param);
Gs = gsp_estimate_lmax(Gs);

% Test with 1d signal
ind = Gs{2}.mr.idx;
fs = rand(size(ind,1),1);

paramint.regularize_epsilon = lambda;
paramint.order = 100;
f1 = zeros(G.N,1);
f1(ind) = fs;
f1 = gsp_interpolate_old(f1,G,ind,paramint);
f2 = gsp_interpolate(G,fs,ind,paramint);

if sum(abs(f1-f2)) < 1e-10
    fprintf('PYRAMID: test interpolate 1d ok\n');
else
    warning('PYRAMID: error in test interpolate 1d');
    badness = sum(abs(f1-f2))
    errors = errors+1;
end


% Test with 2d signal
fs2 = rand(size(ind,1),3);
f3 = gsp_interpolate(G,fs2,ind,paramint);
f4 = zeros(N,3);
for ii = 1:3;
    f4(:, ii) = gsp_interpolate(G,fs2(:, ii),ind,paramint);
end

if f3 == f4;
    fprintf('PYRAMID: test interpolate 2d ok\n');
else
    warning('PYRAMID: error in test interpolate 2d');
    badness = sum(abs(f3-f4))
    errors = errors+1;
end

end



function errors = test_reconstruction()

errors = 0;

N = 256;
Nl = 4;
G = gsp_sensor(N);
Gs = gsp_graph_multiresolution(G,Nl);
Gs = gsp_compute_fourier_basis(Gs);

% Test with 1d signal
f = rand(N,1);
[ca,pe] = gsp_pyramid_analysis(Gs, f,Nl);
% coeff = gsp_pyramid_cell2coeff(ca,pe);
f_pred = gsp_pyramid_synthesis(Gs,ca{end},pe);

if norm(f-f_pred) < eps(100)
    fprintf('PYRAMID: test reconstruction 1d ok\n');
else
    warning('PYRAMID: error in test reconstruction 1d');
    errors = errors+1;
end

% Test with 2d signal
f2 = rand(N,5);
[ca,pe]=gsp_pyramid_analysis(Gs, f2,Nl);
coeff = gsp_pyramid_cell2coeff(ca,pe);
f2_pred = gsp_pyramid_synthesis(Gs,ca{end},pe);

if norm(f2-f2_pred,'fro') < eps(100)
    fprintf('PYRAMID: test reconstruction 2d ok\n');
else
    warning('PYRAMID: error in test reconstruction 2d');
    errors = errors+1;
end

end


function errors = test_result_graph()

errors = 0;

N = 100;
Nl = 1;
G = gsp_sensor(N);
param.sparsify = 0;
gsp_reset_seed
Gs1 = gsp_graph_multiresolution(G,Nl,param);
gsp_reset_seed
Gs2 = gsp_graph_multiresolution_old(G,Nl,param);

if sum(sum(abs(Gs1{2}.L(:)-Gs2{2}.L(:)))) < 1e-12
    fprintf('PYRAMID: test graph 1 ok\n');
else
    warning('PYRAMID: error in test graph 1 ');
    badness = sum(sum(abs(Gs1{2}.L(:)-Gs2{2}.L(:))))
    errors = errors+1;
end


end



function errors = test_result_graph2()

errors = 0;

N = 100;
Nl = 1;
G = gsp_sensor(N);
param.sparsify = 1;
param.sparsify_epsilon = 0.5;
gsp_reset_seed(0);
Gs1 = gsp_graph_multiresolution(G,Nl,param);
gsp_reset_seed(0);
Gs2 = gsp_graph_multiresolution_old(G,Nl,param);

if sum(sum(abs(Gs1{2}.L(:)-Gs2{2}.L(:)))) < 1e-7
    fprintf('PYRAMID: test graph 2 ok\n');
else
    warning('PYRAMID: error in test graph 2 ');
    badness = sum(sum(abs(Gs1{2}.L(:)-Gs2{2}.L(:))))
    errors = errors+1;
end


end


function errors = test_result_graph_multilevel()

errors = 0;

N = 100;
Nl = 3;
G = gsp_sensor(N);
param.sparsify = 0;
param.lambda = 0;
% param.sparsify_epsilon = 0.5;
% param.epsilon = 0.5;
gsp_reset_seed(0)
Gs1 = gsp_graph_multiresolution(G,Nl,param);
gsp_reset_seed(0)
Gs2 = gsp_graph_multiresolution_old(G,Nl,param);

if sum(sum(abs(Gs1{4}.L(:)-Gs2{4}.L(:)))) < 1e-10
    fprintf('PYRAMID: test graph multilevel ok\n');
else
    warning('PYRAMID: error in test graph multilevel ');
    badness = sum(sum(abs(Gs1{4}.L(:)-Gs2{4}.L(:))))
    errors = errors+1;
end


end




function errors = test_result_graph_multilevel2()

errors = 0;

N = 100;
Nl = 3;
G = gsp_sensor(N);
param.sparsify = 1;
param.sparsify_epsilon = 0.5;
param.lambda = 0;
gsp_reset_seed(3)
Gs1 = gsp_graph_multiresolution(G,Nl,param);
gsp_reset_seed(3)
Gs2 = gsp_graph_multiresolution_old(G,Nl,param);
try
if sum(sum(abs(Gs1{4}.L(:)-Gs2{4}.L(:)))) < 1e-7
    fprintf('PYRAMID: test graph multilevel 2 ok\n');
else
    warning('PYRAMID: error in test graph multilevel 2');
    badness = sum(sum(abs(Gs1{3}.L(:)-Gs2{3}.L(:))))
    errors = errors+1;
end
catch
    warning('This test fail but this can happen')
    errors = errors+1;
end

end
