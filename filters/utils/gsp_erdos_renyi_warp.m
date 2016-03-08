function [ cdf ] = gsp_erdos_renyi_warp( ss , N , p )


cutoff=4; % parameter to approximate the non-compactly supported distribution (the free convolution) by a compactly supported one, just for numerical purposes
% To do: pass a parameter
%
%   Url: http://lts2research.epfl.ch/gsp/doc/filters/utils/gsp_erdos_renyi_warp.php

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

num_pts=length(ss);
if num_pts > 2
    delta=min(ss(2:num_pts)-ss(1:(num_pts-1)));
else
    delta=.1;
end

cdf=zeros(size(ss));

for k=1:num_pts
    if ss(k) > (p*N-cutoff*sqrt(p*(1-p)*N)-delta)
        if ss(k) <= (p*N+cutoff*sqrt(p*(1-p)*N)+delta)
            xx=(p*N-cutoff*sqrt(p*(1-p)*N)-2*delta):delta:ss(k);
            cdf(k)=trapz(xx,sqrt(1/((1-p)*N*p))*gsp_free_conv_norm_semi((xx-p*N)/(sqrt(p*(1-p)*N))));
        else
            cdf(k)=1;
        end
    end
end

end


