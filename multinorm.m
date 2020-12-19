function [y return_flag] = multinorm(x,m,covar)
% Evaluates a multidimensional Gaussian
% of mean m and covariance matrix covar
% at the array of points x
% This function ensures that the return value y is valid, so actually the
% return_flag is not used.

real_min = sqrt(realmin);
real_max = sqrt(realmax);

[dim npoints] = size(x);

for i = 1:dim
    if covar(i, i) < real_min
        covar(i, i) =  real_min;
    end
end

[in, cov_mat, det_cov, return_flag] = ve_cov_Jain(covar); % Protection
if return_flag ~= 0 
    y = 0;
    return;
end
ff = ((2*pi)^(-dim/2))*((det_cov)^(-0.5));

quadform = zeros(1,npoints);
centered = (x - m*ones(1,npoints));
if dim ~= 1
   y = ff * exp(-0.5*sum(centered.*(in*centered)));
else
   y = ff * exp(-0.5*in*centered.^2 );
end
return_flag = 0;



