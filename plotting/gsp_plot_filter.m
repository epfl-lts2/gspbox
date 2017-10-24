function [] = gsp_plot_filter(G, filters, param)
%GSP_PLOT_FILTER  Plot a system of filters
%   Usage:  [filter_data,test_sum]=gsp_plot_filter(G,filters);
%           [filter_data,test_sum]=gsp_plot_filter(G,filters,param);
%
%   Input parameters:
%       G       : Graph object (Or lmax)
%       filters : Cell array of filters (or single filter)
%       param   : Optional variable containing additional parameters
%   Output parameters:
%       none
%
%
%   'gsp_plot_filters(G, filter, param)' plots a system of graph spectral
%   filters. 
%
%   Example:::
%
%         Nf = 4;
%         G = gsp_sensor(100);
%         G = gsp_estimate_lmax(G);
%         g = gsp_design_mexican_hat(G, Nf);   
%         gsp_plot_filter(G, g); 
%
%
%   Additional parameters
%   ---------------------
%
%   * *param.line_width* : Width of the filter plots (default 4).
%   * *param.npoints* : Number of point where the filters are evaluated
%     (default 1000).
%   * *param.x_tic* : Distance between x-tick labels.
%   * *param.y_tic* : Distance between y-tick labels (default 0.25).
%   * *param.minor_tick* : To show minor tick marks (default 1).
%   * *param.plot_eigenvalues* : To plot black X marks at all eigenvalues
%     of the graph (You need to compute the Fourier basis to use this
%     option). By default the eigenvalues are plot if they are contained in
%     the Graph.
%   * *param.lambda_highlights* : To plot red X marks at highlight
%     eigenvalues (default 0).
%   * *param.x_width* : Width of X marks for the eigenvalues (default 3).
%   * *param.x_size* : Size of X marks for the eigenvalues (default 8).
%   * *param.show_sum* : To plot an extra line showing the sum of the
%     squared magnitudes of the filters (default 1 if there is multiple
%     filters). 
%   * *param.colors_rgb* : To specify the line colors.
%   * *param.cla* : Clear axis (default 1).
%   * *param.yrange* : To specify a range for the y axis.
%   * *param.verbose* : Verbosity level (1 display the warning - 0 no log)
%     (default 1).
%
%   Demos: gsp_demo
%

% Author : David I Shuman, Nathanael Perraudin
% Testing: test_filter 

  
% Read input parameters
if nargin < 3
   param = struct;
end

if ~isstruct(G)
   G.lmax = G;
   param.plot_eigenvalues = 0;
end

if ~isfield(G,'lmax')
    G = gsp_estimate_lmax(G);
    warning(['GSP_FILTER_ANALYSIS: The variable lmax is not ',...
        'available. The function will compute it for you. ',...
        'However, if you apply many time this function, you ',...
        'should precompute it using the function: ',...
        'gsp_estimate_lmax']);
end

if ~isfield(param,'verbose'), param.verbose = 1; end
if ~isfield(param,'line_width'), param.line_width = 4; end
if ~isfield(param,'x_width'), param.x_width = 3; end
if ~isfield(param,'x_size'), param.x_size = 8; end
if ~isfield(param,'show_sum'), param.show_sum = length(filters)>1; end
if ~isfield(param,'y_tic'), param.y_tic = 0.25; end
if ~isfield(param,'minor_tick'), param.minor_tick = 1; end
if ~isfield(param,'npoints'), param.npoints = 1000; end
if ~isfield(param,'cla'), param.cla = 1; end
if ~isfield(param,'x_tic')
    param.x_tic=max(1,ceil(G.lmax/10));
end
if ~isfield(param,'plot_eigenvalues')
    param.plot_eigenvalues = isfield(G,'e'); 
end

if ~isfield(G,'lmax')
    if param.verbose
        fprintf('GSP_KERNEL_MEXICAN_HAT has to compute lmax \n')
    end
    G = gsp_estimate_lmax(G);
end

lambdas = linspace(0,G.lmax,param.npoints);


if param.cla
    cla;
end
hold on;
% apply the filter
fd = gsp_filter_evaluate(filters,lambdas);

% plot the filter
plot(lambdas,fd,'LineWidth',param.line_width);

% plot the eigenvalues
if param.plot_eigenvalues
    if isfield(G,'e')
        plot(G.e,zeros(G.N,1),'xk','LineWidth',...
            param.x_width,'MarkerSize',param.x_size);
    else
        if param.verbose
            warning('GSP_PLOT_FILTER: No eigenvalues found in the graph');
        end
    end
end

% plot hightlights eigenvalues
if isfield(param,'lambda_highlights')
    plot(param.lambda_highlights, ...
        zeros(length(param.lambda_highlights),1),...
        'xr','LineWidth',param.x_width,'MarkerSize',param.x_size);
end

% plot the sum
if param.show_sum
    test_sum=sum(fd.^2,2); 
    plot(lambdas,test_sum,'k','LineWidth',param.line_width);
end


box on;
% X axis
xlim(full([0 G.lmax]));
% Y axis
if isfield(param,'yrange')
    yrange=param.yrange;
    ylim(yrange);
    set(gca,'YTick',yrange(1):y_tic:yrange(2));
end

% Add marks
if param.minor_tick
    set(gca,'XMinorTick','on','YMinorTick','on');

end

% Change the color
if isfield(param,'colors_rgb');
    set(gca, 'ColorOrder', param.colors_rgb);
end
set(gca,'XTick',0:param.x_tic:G.lmax);


hold off;
end


