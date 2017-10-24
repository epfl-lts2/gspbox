function out = gsp_filter(G, fi, signal, param)
%GSP_FILTER Filter function
%   Usage:  coeffs = gsp_filter(G, fi, signal);
%           coeffs = gsp_filter(G, fi, signal, param);
%
%   Input parameters:
%         G         : Graph structure.
%         fi        : Spectral filter.
%         s         : Graph signal to filters
%         param     : Optional parameters
%   Output parameters:
%         c         : Filtered signal
%
%   This function is a shortcut to the function |gsp_filter_analysis|.
%   Please use the documentation of |gsp_filter_analysis|
%   

if nargin<4
    param = struct;
end

    out = gsp_filter_analysis(G, fi, signal, param);
end
