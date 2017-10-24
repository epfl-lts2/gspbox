function [interpolated_values] = gsp_pwl_warp_fn(x,y,x0)

cut=1e-4;

if max(x0)>max(x)+cut
    error('GSP_PWL_WARP_FN: This function does not allow you to interpolate ousite the point x and y')
end

if min(x0)<min(x)-cut
    error('GSP_PWL_WARP_FN: This function does not allow you to interpolate ousite the point x and y')
end

sx0 = size(x0);
x0 = x0(:);


% Make sure data is sorted and monotonic
[x,x_ind]=sort(x,'ascend');
y=y(x_ind);
if ( isequal(sort(y,'ascend'),y)==0 && isequal(sort(y,'descend'),y)==0 )
    error('Data points are not monotonic');
end

% Piecewise-linear interpolation
num_pts=length(x);
num_pts_to_interpolate=length(x0);
interpolated_values=zeros(num_pts_to_interpolate,1);

for i=1:num_pts_to_interpolate
    [~,closest_ind]=min(abs(x-x0(i)));
    if (x(closest_ind)-x0(i))<(-cut) || ( abs(x(closest_ind)-x0(i))<cut && closest_ind < num_pts)
    %if sign(x(closest_ind)-x0(i))<0 || ( sign(x(closest_ind)-x0(i))==0 && closest_ind < num_pts)
        lower_ind=closest_ind;
    else
        lower_ind=closest_ind-1;
    end
    interpolated_values(i)=y(lower_ind)*(x(lower_ind+1)-x0(i))/(x(lower_ind+1)-x(lower_ind)) + y(lower_ind+1)*(x0(i)-x(lower_ind))/(x(lower_ind+1)-x(lower_ind));
end

interpolated_values = reshape(interpolated_values,sx0);

end

% Old Code:
% function [output] = gsp_pwl_warp_fn(approx_spectrum,x0)
% N=length(approx_spectrum);
% output=zeros(size(x0));
% for i=1:length(x0);
%     left_index=sum(approx_spectrum<=x0(i));
%     left_val=approx_spectrum(max(left_index,1));
%     right_val=approx_spectrum(min(left_index+1,N));
%     output(i)=left_index+(x0(i)-left_val)/(right_val-left_val);
% end
% 
% end

