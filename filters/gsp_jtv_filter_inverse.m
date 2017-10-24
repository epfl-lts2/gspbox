function [s] = gsp_jtv_filter_inverse(G, filter, filtertype,c, param)
%GSP_JTV_FILTER_INVERSE Inverse operator of a joint timve_vertex filterbank
%   Usage:  s = gsp_jtv_filter_inverse(G, filter, c);
%           s = gsp_jtv_filter_inverse(G, filter, c, param);
%
%   Input parameters:
%         G          : Time-vertex Graph structure.
%         filter     : Cell array of time-vertex filters.
%         filtertype : Filter domain (ts,js,ts-array,js-array)
%         c          : Transform coefficients
%         param      : Optional parameter
%   Output parameters:
%         signal     : sythesis signal
%

% Author: Francesco Grassi
% Testing: test_jtv_filter
% Date: September 2016

if nargin<5
    param = struct;
end

[dual_filter,filtertype] = gsp_jtv_design_can_dual(filter,filtertype);

s = gsp_jtv_filter_synthesis(G,dual_filter,filtertype,c,param);


end



