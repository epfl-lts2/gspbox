function fd = gsp_filter_evaluate(filter, x)
%GSP_FILTER_EVALUATE Evaluate the filterbank
%   Usage: fd = gsp_filter_evaluate(filter, x)
%
%   Input parameters:
%       filter  : cell array of filter
%       x       : data
%   Output parameters:
%       fd      : response of the filters
%
%   This function applies all the filters in *filter* to the data *x*. Every
%   filter correspond to one column of the matrix *fd*.
%
%   See also: gsp_filter_analysis
%

% Author: Nathanael Perraudin
% Date: 18 March 2014

if ~iscell(filter)
   filter = {filter};
end

Nf=numel(filter);
fd=zeros(length(x),Nf);
for k=1:Nf
    fd(:,k)=filter{k}(x);
end

end
