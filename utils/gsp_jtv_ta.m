function [t,label] = gsp_jtv_ta(G,lag)
%GSP_JTV_TA Time axis for the joint time vertex framework
%   Usage:  t = gsp_jtv_ta(G);
%
%   Input parameters:
%       G       : Time-vertex graph structure
%       lag     : Compute lag axis
%   Ouput parameters:
%       t       : Time or Lag axis (row vector)

% Author: Francesco Grassi
% Date  : September 2016

if nargin<2
    lag = 0;
end

if isfield(G.jtv,'T')
    T = G.jtv.T;
    fs = G.jtv.fs;
    if lag
        t = -(T-1)/fs:1/fs:(T-1)/fs;
        label = 'Lag';
    else
        t = 0:1/fs:(T-1)/fs;
        label = 'Time';
    end
else
    error('GSP_JTV_TA needs the time dimension. Use GSP_JTV_GRAPH');
end

end