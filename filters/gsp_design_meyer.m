function [ g,t ] = gsp_design_meyer(G, Nf, param)
%GSP_DESIGN_MEYER Design the meyer filterbank
%   Usage: g =  gsp_design_meyer(G, Nf, param);
%               gsp_design_meyer(G ,Nf);
%               gsp_design_meyer(G);
%
%   Input parameters:
%         G             : Graph or upper bound on the Laplacian spectrum
%         Nf            : Number of filters to cover the interval [0,lmax] (default 6)
%         param         : Structure of optional parameters
%   Output parameters:
%         g             : A cell array of filters
%
%   This function return a array of filters designed to be meyer wavelet.
%
%   *param* is an optional structure containing the following fields
%
%   * *param.t*: vector of scale to be used (default: log scale)
%   * *param.verbose*: verbosity level. 0 no log - 1 display warnings.
%     (default 1) 
%
%   This function will compute the maximum eigenvalue of the laplacian. To
%   be more efficient, you can precompute it using::
%
%       G = gsp_estimate_lmax(G);
%
%   Example:::
%
%         Nf = 4;
%         G = gsp_sensor(100);
%         G = gsp_estimate_lmax(G);
%         g = gsp_design_meyer(G, Nf);   
%         gsp_plot_filter(G,g);  
%
%   This function is inspired by the sgwt_toolbox. 
%       
%   See also:

% Author: Nathanael Perraudin, David K. Hammond
% Date: 18 March 2014

    

if nargin < 3
    param = struct;
end

if ~isfield(param,'verbose'), param.verbose = 1; end

if isstruct(G)
    if ~isfield(G,'lmax')
        if param.verbose
            fprintf('GSP_DESIGN_MEYER has to compute lmax \n')
        end
        G = gsp_estimate_lmax(G);
    end
   lmax = G.lmax;
else
   lmax = G;
end

if ~isfield(param,'t')
   param.t = (4/(3*lmax)) * 2.^(Nf-2:-1:0);
end


if param.verbose
    if length(param.t) ~= Nf - 1
       warning(['GSP_KERNEL_MEYER: You have specified ',...
           'more scales than Number of filter -1']);
    end
end

t = param.t;
g = cell(Nf,1);

g{1}= @(x) kernel_meyer(t(1)*x,'sf');
for j=1:Nf-1
    g{j+1}= @(x) kernel_meyer(t(j)*x,'wavelet');
end



end




function r=kernel_meyer(x,kerneltype)
% sgwt_kernel_meyer : evaluates meyer wavelet kernel and scaling function
% function r=sgwt_kernel_meyer(x,kerneltype)
%
% Inputs
% x : array of independent variable values
% kerneltype : string, either 'sf' or 'wavelet' 
%
% Ouputs
% r : array of function values, same size as x.
%
% meyer wavelet kernel : supported on [2/3,8/3]
% meyer scaling function kernel : supported on [0,4/3]
%
% Use of this kernel for SGWT proposed by Nora Leonardi and Dimitri Van De Ville,
% "Wavelet Frames on Graphs Defined by fMRI Functional Connectivity"
% International Symposium on Biomedical Imaging, 2011
l1=2/3;
l2=4/3;%2*l1;
l3=8/3;%4*l1;
v=@(x) x.^4.*(35-84*x+70*x.^2-20*x.^3) ; 

% as we initialize r with zero, computed function will implicitly be zero for
% all x not in one of the three regions defined above
r=zeros(size(x));
switch kerneltype
  case 'sf'
    r(x<l1)=1;
    r(x>=l1 & x<l2)=cos((pi/2)*v(abs(x(x>=l1 & x<l2))/l1-1));
  case 'wavelet'
    r(x>=l1 & x<l2)=sin((pi/2)*v(abs(x(x>=l1 & x<l2))/l1-1));
    r(x>=l2 & x<l3)=cos((pi/2)*v(abs(x(x>=l2 & x<l3))/l2-1));
  otherwise
    error(sprintf('unknown kernel type %s',kerneltype));
end


end

