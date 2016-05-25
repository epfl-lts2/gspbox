function [ interpolated_values ] = gsp_mono_cubic_warp_fn(x,y,x0)

cut=1e-4;

% Make sure data is sorted and monotonic
%
%   Url: http://lts2research.epfl.ch/gsp/doc/filters/utils/gsp_mono_cubic_warp_fn.php

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.6.0
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
[x,x_ind]=sort(x,'ascend');
y=y(x_ind);
if ( isequal(sort(y,'ascend'),y)==0 && isequal(sort(y,'descend'),y)==0 )
    error('Data points are not monotonic');
end

% Monotonic cubic interpolation using the Fritsch-Carlson method
num_pts=length(x);
if length(y) ~= num_pts
    error('x and y vectors have different dimensions');
end

% 1. Compute slopes of secant lines
Delta=(y(2:end)-y(1:num_pts-1))./(x(2:end)-x(1:num_pts-1));

% 2. Initialize tangents m at every data point
m = (Delta(1:num_pts-2)+Delta(2:num_pts-1))/2;
m = [Delta(1);m;Delta(end)];

% 3. Check for equal y's to set slopes equal to zero
for k=1:num_pts-1
    if Delta(k)==0
        m(k)=0;
        m(k+1)=0;
    end
end

% 4. Initialize alpha and beta
alpha = m(1:num_pts-1)./Delta;
beta = m(2:num_pts)./Delta;

% 5. Make monotonic
for k=1:num_pts-1
    if alpha(k)^2+beta(k)^2 > 9
        tau=3/sqrt(alpha(k)^2+beta(k)^2);
        m(k)=tau*alpha(k)*Delta(k);
        m(k+1)=tau*beta(k)*Delta(k);
    end
end

% 6. Cubic interpolation
num_pts_to_interpolate=length(x0);
interpolated_values=zeros(size(x0));

for ii=1:num_pts_to_interpolate
    [~,closest_ind]=min(abs(x-x0(ii)));
    %if sign(x(closest_ind)-x0(i))<0 || ( sign(x(closest_ind)-x0(i))==0 && closest_ind < num_pts)
    if (x(closest_ind)-x0(ii))<-cut || ( abs(x(closest_ind)-x0(ii))<cut && closest_ind < num_pts)
        lower_ind=closest_ind;
    else
        lower_ind=closest_ind-1;
    end
    h=x(lower_ind+1)-x(lower_ind);
    t=(x0(ii)-x(lower_ind))/h;

    interpolated_values(ii) = y(lower_ind)*(2*t^3-3*t^2+1) + h*m(lower_ind)*(t^3-2*t^2+t) + y(lower_ind+1)*(-2*t^3+3*t^2) + h*m(lower_ind+1)*(t^3-t^2); 
end
end


