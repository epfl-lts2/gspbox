function [G, pixels, patches] = gsp_patch_graph(img, param)
%GSP_PATCH_GRAPH Create a graph by NN patches of an image
%   Usage :  G = compute_patch_graph( img );
%            [G, f] = gsp_patch_graph( img, param );
%            [G, f, patches] = gsp_patch_graph( img, param );
%
%   Input parameters:
%       img         : Input image (unknown pixels marked as NaN)
%       param       : Structure of optional parameters
%
%   Output parameters:
%       G           : Resulting graph
%       f           : Image signal
%       patches     : Patches
%
%   'compute_patch_graph( path , param )' creates a graph between pixels in
%   an image by connecting them using the euclidean distance between
%   patches around the pixels
%
%   Additional parameters
%   ---------------------
%
%   * *param.patch_size*      : int     the patch size in pixels (odd)
%   * *param.scale*           : float   to rescale the input image
%   * *param.rho*             : float   spatial constraint
%   * *use_incomplete_patch*  : boolean use incomplete patch for the graph construction
%   * *param.nnparam*         : struct  parameters for graph construction
%

% Author: Johan Paratte, Michael Defferrard, Nathanael Perraudin
% Date: November 2014
% 

    if nargin < 2
    % Define parameters
        param = {};
    end
    
    if ~isfield(param, 'patch_size'), param.patch_size = 5; end
    if ~isfield(param, 'rho'), param.rho = 0.001; end
    if ~isfield(param, 'use_incomplete_patch'), param.use_incomplete_patch = 0; end
    if ~isfield(param, 'nnparam')
       param.nnparam = {};
       param.nnparam.center = 0;
       param.nnparam.resize = 0;
       param.nnparam.k = 10;
       param.nnparam.use_l1 = 0;
       param.nnparam.use_flann = 1;
    end
    
    if (mod(param.patch_size, 2) == 0) 
        disp('Patch size must be odd, converting to closest odd number');
        param.patch_size = param.patch_size + 1;
    end
    
%     if length(size(img)) > 2
%        img = rgb2gray(img); 
%     end

    %Extract patches
    [oheight, owidth,Nc] = size(img);
    
    
    
    psize = param.patch_size;
    
    %Parameters
    if isfield(param, 'scale')
        img = imresize(img, param.scale);
    end
    
    margin = floor(psize / 2);
    

    
    dim = psize*psize;
    %height = (oheight - 2*margin);
    %width = (owidth - 2*margin);
    
    %Expand the image to compensate the margin
    h = oheight;
    w = owidth;
    new_img = zeros(h + 2*margin, w + 2*margin,Nc);
    %Copy image
    new_img(1+margin:end-margin, 1+margin:end-margin,:) = img;
    %Four lines
    new_img(1:margin, 1+margin:end-margin,:) = img(margin:-1:1, :,:);
    new_img(1+margin:end-margin, 1:margin,:) = img(:, margin:-1:1,:);
    new_img(end-margin:end, 1+margin:end-margin,:) = img(end:-1:end-margin, :,:);
    new_img(1+margin:end-margin, end-margin:end,:) = img(:,end:-1:end-margin,:);
    %Four corners
    new_img(1:margin, 1:margin,:) = img(margin:-1:1, margin:-1:1,:);
    new_img(1:margin, end-margin:end,:) = img(margin:-1:1, end:-1:end-margin,:);
    new_img(end-margin:end, 1:margin,:) = img(end:-1:end-margin, margin:-1:1,:);
    new_img(end-margin:end, end:-1:end-margin,:) = img(end:-1:end-margin, end:-1:end-margin,:);
    
    % Signals on the graph.
    patches = zeros(h*w, dim*Nc+2);
    pixels  = double(reshape(img, h*w, Nc)); % zeros(h*w, 1);
    coords  = zeros(h*w, 2);  % gsp_plot_signal takes [x,y]
    coords(:,2) = repmat((h:-1:1).', w, 1);
    coords(:,1) = reshape(repmat((1:+1:w), h, 1), [], 1);
    
    % Pre-allocation.
    nUnknown = sum(sum(img(:,:,1)<0)) * 2;
    unknowns = zeros(nUnknown, 1);
    iUnknown = 0;
    
    count = 1;
    
    % For speed improvement, one of the for loop can probably be removed
    % here
    for w = 1:owidth
       for h = 1:oheight
            patches(count, 1:Nc*dim) = reshape(new_img(h:h+psize-1, w:w+psize-1,:), 1, Nc*dim);
            patches(count, Nc*dim+1:end) = [param.rho*w; param.rho*h];
            count = count + 1;
        end
    end
    
    fprintf('Compute graph with %d vertices\n', size(patches,1));
    
    if param.use_incomplete_patch

        G = gsp_rmse_mv_graph(patches, param.nnparam);
    else
        % Ignore patches which contain unknown pixel values.
        % We do not want them to be connected with known patches.
        clearedPatches = patches;

        for count = 1:size(clearedPatches,1)
            if sum(isnan(clearedPatches(count,:)))
                clearedPatches(count,:) = -1e3;
                % List of unknown patches.
                iUnknown = iUnknown + 1;
                unknowns(iUnknown) = count;
            end
        end


        G = gsp_nn_graph(clearedPatches, param.nnparam);

        % Disconnect unknown patches.
        unknowns = unknowns(1:iUnknown);
        G.W(unknowns,:) = 0;
        G.W(:,unknowns) = 0;
        
        % Update G.A, G.d and G.Ne.
        G = gsp_graph_default_parameters(G);
    end
    

    
    % These are coordinates in the high dimensional patch space. We prefer
    % to visualize the graph in the 2D image space.
%     patches = G.coords;
    G.coords = coords;
    G.plotting.limits = [1,w,1,h];
    
end
