function [ cdf ] = gsp_erdos_renyi_warp( ss , N , p )


cutoff=4; % parameter to approximate the non-compactly supported distribution (the free convolution) by a compactly supported one, just for numerical purposes
% To do: pass a parameter

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

