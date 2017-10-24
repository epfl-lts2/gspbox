%GSP_DEMO Tutorial on the GSPBox
% 
%   In this demo, we are going to show the basic operations of the GSPBox.
%   To lauch the toolbox, just go into the repository where the GSPBox was
%   extracted and type:
%
%           gsp_start;
%
%   A banner will popup telling you that everything happens correctly. To
%   speedup some processing, you might want to compile some mexfile. Refer
%   to |gsp_make| for more informations. However, if the compilation is not
%   working on your computer, keep quiet, everything should still work and
%   most of the routine are implemented only in matlab.
%
%   Most likely, the first thing you would like to do is to create a graph.
%   To do so, you only need the adjacendy or the weight matrix $W$. Once
%   you have it, you can construct a graph using::
%
%           G = gsp_graph(W);
%
%   This function will create a full structure ready to be used with the
%   toolbox. To know a bit more about what is in this structure, you can
%   refer to the help of the function |gsp_graph_default_parameters|.
%
%   The GSPBox contains also a list of graph generators. To see a full list
%   of these graphs, type:::
%
%           help graphs
%
%   For this demo, we will use the graph |gsp_logo|. You can load it
%   using:::
%
%           G = gsp_logo
%
%   Here observe the attribute of the structure *G*. 
%
%   * *G.W*: Weight matrix 
%   * *G.A*: Adacency matrix 
%   * *G.N*: Number of nodes 
%   * *G.type*: Type of graph 
%   * *G.directed*: 1 if the graph is directed, 0 if not
%   * *G.lap_type*: Laplacian type 
%   * *G.d*: Degree vector 
%   * *G.Ne*: Number of edges
%   * *G.coords*: Coordinates of the vertices
%   * *G.plotting*: Plotting parameters 
%
%   In the folder 'plotting', the GSPBox contains some plotting routine.
%   For instance, we can plot a graph using::
%
%           gsp_plot_graph(G);
%
%   .. figure::
%
%      GSP graph
%
%      This figure shows the result of the command 'gsp_plot_graph(G)'
%
%   Wonderful! Isn't it? Now, let us start to analyse this graph. To compute
%   graph Fourier transform or exact graph filtering, you need to
%   precompute the Fourier basis of the graph. This operation could be
%   relatively long since it involves a full diagonalization of the
%   Laplacian. Don't worry, you do not need to perform this operation to
%   filter signals on graph. The fourier basis is computed by::
%
%           G = gsp_compute_fourier_basis(G);
%
%   The function |gsp_compute_fourier_basis| add two new fields to the
%   structure *G*:
%
%   * *G.U*: The eigenvectors of the Fourier basis
%   * *G.e*: The eigenvalues
%
%   The fourier eigenvectors does look like a sinusoide on the graph. Let's
%   plot the second and the third ones. (The first one is constant!)::
%
%           gsp_plot_signal(G,G.U(:,2));
%           title('Second eigenvector')
%           subplot(212)
%           gsp_plot_signal(G,G.U(:,3));
%           title('Third eigenvector')
%
%   .. figure::
%
%      Eigenvectors
%
%
%
%   Now, we are going to show a basic filtering operation. Filters are usually
%   defined in the spectral domain. To define the following filter
%
%   ..   h(x) = 1/(1+tau*x),
%
%   .. math:: h(x) =\frac{1}{1+\tau x},
%
%   just write in Matlab::
%
%           tau = 1;
%           h = @(x) 1./(1+tau*x);
%
%   Hint: You can define filterbank using cell array!
%
%   Let's display this filter::
%
%           gsp_plot_filter(G,h);
%
%   .. figure::
%
%      Low pass filter $h$
%
%      The filter $h$ is plotted along all the spectrum of the graph.
%      The black cross are the eigenvalues of the Laplacian. They are the
%      points where the continuous filter will be evaluated to create a
%      discrete filter.
%
%   To apply the filter to a given signal, you only need to run a single
%   function::
%
%           % Create a signal
%           f = zeros(G.N,1);
%           f(G.info.idx_g) = -1;
%           f(G.info.idx_s) = 1;
%           f(G.info.idx_p) = -0.5;
%           f = f + 0.3*randn(G.N,1);
%           % Remove the noise
%           f2 = gsp_filter(G,h,f);
%
%   `gsp_filter` is actually a shortcut to |gsp_filter_analysis|.
%   `gsp_filter_analysis` performs the analysis operator associated to a
%   filterbank. See the |gsp_demo_wavelet| for more information.
%
%   Finnaly, we display the result of this low pass filtering on the graph::
%
%           figure;
%           subplot(211)
%           gsp_plot_signal(G,f);
%           title('Signal with noise')
%           subplot(212)
%           gsp_plot_signal(G,f2);
%           title('Signal denoised');
%
%   .. figure::
%
%      Result of filtering
%
%      The noise is largely removed thanks to the filter. However, some
%      energy is diffused between the letters. This is the typical
%      behaviour of a low pass filter.
%
%   Enjoy the GSPBOX !
%


% Author: Nathanael Perraudin
% Date : 14 August 2014

clear;
close all;

G = gsp_logo;

% display the graph
figure;
gsp_plot_graph(G);

%% Compute the Fourier basis

G = gsp_compute_fourier_basis(G);

% Display an eigenvector
figure;
subplot(211)
gsp_plot_signal(G,G.U(:,2));
title('Second eigenvector')
subplot(212)
gsp_plot_signal(G,G.U(:,3));
title('Third eigenvector')

%% Create a signal

f = zeros(G.N,1);
f(G.info.idx_g) = -1;
f(G.info.idx_s) = 1;
f(G.info.idx_p) = -0.5;
f = f + 0.3*randn(G.N,1);


%% Define a low pass filter
tau = 1;
h = @(x) 1./(1+tau*x);

figure
gsp_plot_filter(G,h);
title('Filter h')


%% Perform the filtering operation
f2 = gsp_filter(G,h,f);
% f2 = gsp_filter_analysis(G,h,f);


% Display the result
figure;
subplot(211)
gsp_plot_signal(G,f);
title('Signal with noise')
subplot(212)
gsp_plot_signal(G,f2);
title('Signal denoised');



