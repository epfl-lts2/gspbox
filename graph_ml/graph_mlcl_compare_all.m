function [rel_error, name ] = graph_mlcl_compare_all(x,xx,y,yy, param)
%GRAPH_ML_COMPARE_ALL Compare all methods for a classification problem
%
%   Url: http://lts2research.epfl.ch/gsp/doc/graph_ml/graph_mlcl_compare_all.php

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.5.1
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% If you use this toolbox please kindly cite
%     N. Perraudin, J. Paratte, D. Shuman, V. Kalofolias, P. Vandergheynst,
%     and D. K. Hammond. GSPBOX: A toolbox for signal processing on graphs.
%     ArXiv e-prints, Aug. 2014.
% http://arxiv.org/abs/1408.5781

if nargin<5
    param = struct;
end

if ~isfield(param,'verbose'), param.verbose = 1; end
if ~isfield(param,'tol'), param.tol = 1e-8; end
if ~isfield(param,'k'), param.k = 10; end
if ~isfield(param,'maxit'), param.maxit = 2000; end


%% Error function
Nl = size(x, 1);
ev = @(x) sum(sum(abs(x((Nl+1):end)-yy)>0))/numel(yy);
rel_error = [];
name = {};
%% Prepare data

xtot = [x; xx];


M = zeros(size(xtot,1),1);
M(1:size(x, 1)) = 1;
ytot = [y; yy];
ytot = M.*ytot;

if param.verbose
    fprintf('Start the simulations \n ');
end







%% KNN classifier
methodname = 'KNN Classifier';

clear paramm

paramm = param;
paramm.use_flann = 0;
paramm.weighted = 0;

Gcknn = gsp_knn_classify_graph(x,xx,paramm);
s = gsp_classification_knn(Gcknn,M,ytot);

[rel_error, name ] = update_values(ev(s), rel_error, name , methodname,param);

%% KNN classifier FLANN
methodname = 'KNN Classifier FLANN';

clear paramm

paramm = param;
paramm.use_flann = 1;
paramm.weighted = 1;

Gcknnflann = gsp_knn_classify_graph(x,xx,paramm);
s = gsp_classification_knn(Gcknnflann,M,ytot);

[rel_error, name ] = update_values(ev(s), rel_error, name , methodname, param);


%% KNN classifier weighted
methodname = 'KNN Classifier weighted';

clear paramm

paramm = param;
paramm.use_flann = 0;
paramm.weighted = 1;


Gcknnw = gsp_knn_classify_graph(x,xx,paramm);
s = gsp_classification_knn(Gcknnw,M,ytot);

[rel_error, name ] = update_values(ev(s), rel_error, name , methodname,param);

%% KNN classifier weighted FLANN
methodname = 'KNN Classifier weighted FLANN';

clear paramm

paramm = param;
paramm.use_flann = 1;
paramm.weighted = 1;

Gcknnflannw = gsp_knn_classify_graph(x,xx,paramm);
s = gsp_classification_knn(Gcknnflannw,M,ytot);

[rel_error, name ] = update_values(ev(s), rel_error, name , methodname, param);



%% TIK
methodname = 'TIK';

clear paramm

paramm = param;
paramm.use_flann = 0;

G = gsp_nn_graph(xtot,paramm);

s = gsp_classification_tik(G,M,ytot, 0 , param);

[rel_error, name ] = update_values(ev(s), rel_error, name , methodname,param);


%% TIK FLANN
methodname = 'TIK FLANN';

clear paramm

paramm = param;
paramm.use_flann = 1;

Gflann = gsp_nn_graph(xtot,paramm);

s = gsp_classification_tik(Gflann,M,ytot, 0 , param);

[rel_error, name ] = update_values(ev(s), rel_error, name , methodname,param);


%% TIK non weighted
methodname = 'TIK non weighted';

clear paramm

paramm = param;
paramm.use_flann = 0;

Gknn = gsp_nn_graph(xtot,paramm);
Gknn.W = Gknn.W > 0;
Gknn = gsp_graph_default_parameters(Gknn);

s = gsp_classification_tik(Gknn,M,ytot, 0 , param);

[rel_error, name ] = update_values(ev(s), rel_error, name , methodname,param);


%% TIK non weighted FLANN
methodname = 'TIK non weighted FLANN';

clear paramm

paramm = param;
paramm.use_flann = 1;

Gknnflann = gsp_nn_graph(xtot,paramm);
Gknnflann.W = Gknnflann.W > 0;
Gknnflann = gsp_graph_default_parameters(Gknnflann);

s = gsp_classification_tik(Gknnflann,M,ytot, 0 , param);

[rel_error, name ] = update_values(ev(s), rel_error, name , methodname,param);




%% TIK normalized
methodname = 'TIK normalized';

G = gsp_create_laplacian(G,'normalized');
s = gsp_classification_tik(G,M,ytot, 0 , param);

[rel_error, name ] = update_values(ev(s), rel_error, name , methodname,param);


%% TIK normalized FLANN
methodname = 'TIK normalized FLANN';

Gflann = gsp_create_laplacian(Gflann,'normalized');
s = gsp_classification_tik(Gflann,M,ytot, 0 , param);

[rel_error, name ] = update_values(ev(s), rel_error, name , methodname,param);


%% TIK non weighted - normalized
methodname = 'TIK non weighted - normalized';


Gknn = gsp_create_laplacian(Gknn,'normalized');
s = gsp_classification_tik(Gknn,M,ytot, 0 , param);

[rel_error, name ] = update_values(ev(s), rel_error, name , methodname,param);


%% TIK non weighted - normalized FLANN
methodname = 'TIK non weighted - normalized FLANN';

Gknnflann = gsp_create_laplacian(Gknnflann,'normalized');
s = gsp_classification_tik(Gknnflann,M,ytot, 0 , param);

[rel_error, name ] = update_values(ev(s), rel_error, name , methodname,param);

%% TIK tau 0.01
methodname = 'TIK - tau 0.01';

G = gsp_create_laplacian(G,'combinatorial');
s = gsp_classification_tik(G,M,ytot, 0.01*G.Ne/G.N , param);

[rel_error, name ] = update_values(ev(s), rel_error, name , methodname,param);

%% TIK - tau 0.01 normalized
methodname = 'TIK - tau 0.01 normalized';

G = gsp_create_laplacian(G,'normalized');
s = gsp_classification_tik(G,M,ytot, 0.01*G.Ne/G.N  , param);

[rel_error, name ] = update_values(ev(s), rel_error, name , methodname,param);

%% TIK tau 0.1
methodname = 'TIK - tau 0.1';

G = gsp_create_laplacian(G,'combinatorial');
s = gsp_classification_tik(G,M,ytot, 0.1*G.Ne/G.N , param);

[rel_error, name ] = update_values(ev(s), rel_error, name , methodname,param);

%% TIK - tau 0.1 normalized
methodname = 'TIK - tau 0.1 normalized';

G = gsp_create_laplacian(G,'normalized');
s = gsp_classification_tik(G,M,ytot, 0.1*G.Ne/G.N  , param);

[rel_error, name ] = update_values(ev(s), rel_error, name , methodname,param);

%% TIK tau 1
methodname = 'TIK - tau 1';

G = gsp_create_laplacian(G,'combinatorial');
s = gsp_classification_tik(G,M,ytot, 1*G.Ne/G.N , param);

[rel_error, name ] = update_values(ev(s), rel_error, name , methodname,param);

%% TIK - tau 1 normalized
methodname = 'TIK - tau 1 normalized';

G = gsp_create_laplacian(G,'normalized');
s = gsp_classification_tik(G,M,ytot, G.Ne/G.N  , param);

[rel_error, name ] = update_values(ev(s), rel_error, name , methodname,param);

%% TV tau 0.01
methodname = 'TV - tau 0.01';

G = gsp_create_laplacian(G,'combinatorial');
s = gsp_classification_tv(G,M,ytot, 0.01*G.Ne/G.N , param);

[rel_error, name ] = update_values(ev(s), rel_error, name , methodname,param);

%% TV - tau 0.01 normalized
methodname = 'TV - tau 0.01 normalized';

G = gsp_create_laplacian(G,'normalized');
s = gsp_classification_tv(G,M,ytot, 0.01*G.Ne/G.N  , param);

[rel_error, name ] = update_values(ev(s), rel_error, name , methodname,param);

%% TV tau 0.1
methodname = 'TV - tau 0.1';

G = gsp_create_laplacian(G,'combinatorial');
s = gsp_classification_tv(G,M,ytot, 0.1*G.Ne/G.N , param);

[rel_error, name ] = update_values(ev(s), rel_error, name , methodname,param);

%% TV - tau 0.1 normalized
methodname = 'TV - tau 0.1 normalized';

G = gsp_create_laplacian(G,'normalized');
s = gsp_classification_tv(G,M,ytot, 0.1*G.Ne/G.N  , param);

[rel_error, name ] = update_values(ev(s), rel_error, name , methodname,param);

%% TV tau 1
methodname = 'TV - tau 1';

G = gsp_create_laplacian(G,'combinatorial');
s = gsp_classification_tv(G,M,ytot, 1*G.Ne/G.N , param);

[rel_error, name ] = update_values(ev(s), rel_error, name , methodname,param);

%% TV - tau 1 normalized
methodname = 'TV - tau 1 normalized';

G = gsp_create_laplacian(G,'normalized');
s = gsp_classification_tv(G,M,ytot, G.Ne/G.N  , param);

[rel_error, name ] = update_values(ev(s), rel_error, name , methodname,param);


%% TV
methodname = 'TV';

G = gsp_create_laplacian(G,'combinatorial');
s = gsp_classification_tv(G,M,ytot, 0 , param);

[rel_error, name ] = update_values(ev(s), rel_error, name , methodname,param);


%% TV FLANN
methodname = 'TV FLANN';

Gflann = gsp_create_laplacian(Gflann,'combinatorial');
s = gsp_classification_tv(Gflann,M,ytot, 0 , param);

[rel_error, name ] = update_values(ev(s), rel_error, name , methodname,param);


%% TV non weighted
methodname = 'TV non weighted ';


Gknn = gsp_create_laplacian(Gknn,'combinatorial');
s = gsp_classification_tv(Gknn,M,ytot, 0 , param);

[rel_error, name ] = update_values(ev(s), rel_error, name , methodname,param);


%% TV non weighted FLANN
methodname = 'TV non weighted FLANN';

Gknnflann = gsp_create_laplacian(Gknnflann,'combinatorial');
s = gsp_classification_tv(Gknnflann,M,ytot, 0 , param);

[rel_error, name ] = update_values(ev(s), rel_error, name , methodname,param);



%% TV normalized
methodname = 'TV normalized';

G = gsp_create_laplacian(G,'normalized');
s = gsp_classification_tv(G,M,ytot, 0 , param);

[rel_error, name ] = update_values(ev(s), rel_error, name , methodname,param);


%% TV normalized FLANN
methodname = 'TV normalized FLANN';

Gflann = gsp_create_laplacian(Gflann,'normalized');
s = gsp_classification_tv(Gflann,M,ytot, 0 , param);

[rel_error, name ] = update_values(ev(s), rel_error, name , methodname,param);


%% TV non weighted - normalized
methodname = 'TV non weighted - normalized';


Gknn = gsp_create_laplacian(Gknn,'normalized');
s = gsp_classification_tv(Gknn,M,ytot, 0 , param);

[rel_error, name ] = update_values(ev(s), rel_error, name , methodname,param);


%% TV non weighted - normalized FLANN
methodname = 'TV non weighted - normalized FLANN';

Gknnflann = gsp_create_laplacian(Gknnflann,'normalized');
s = gsp_classification_tv(Gknnflann,M,ytot, 0 , param);

[rel_error, name ] = update_values(ev(s), rel_error, name , methodname,param);

%%
if param.verbose
    fprintf('Ends of simulations \n ');
    bar(rel_error)
    set(gca,'XTickLabel',name)

end

%%
end
 



function [rel_error, name ] = update_values(e, rel_error, name , methodname, param)

rel_error = [rel_error; e];
N = length(rel_error);
if N==1
    name = {methodname};
else
    name = {name{1:(N-1)},methodname};
end

if param.verbose
    fprintf([' Sim: %2i,  %3.2f percent -- ', methodname,'\n ' ],N,e*100);
end

end
