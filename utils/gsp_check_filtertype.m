function bool = gsp_check_filtertype(filtertype,type)
%GSP_CHECK_FILTERTYPE Check the type of filter used
%
%   Usage:  bool = gsp_check_fourier(G):
%
%   Input parameters:
%      filtertype           : filtertype to check
%      type                 : type list
%   Output parameters:
%       bool        : boolean
%
%   Check the type of filter used between js, js-array, ts, ts-array
%

if nargin<2
    type = {'js','js-array','ts','ts-array'};
end


bool = 0;

if ~ischar(filtertype)
    return
end

if ~any(strcmpi(filtertype,type))
    return
end

bool = 1;