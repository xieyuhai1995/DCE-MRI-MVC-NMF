% VE_COV

function [inv_cov, cov_mat, det_cov, return_flag] = ve_cov_Jain(cov);
% Actually, there two problems in this ve_cov function
% The first one is the singularity problem; the second one is the realmin
% problem. The singularity problem can be solved by the energy assignment;
% the realmin problem can only be solved by multiply the input data by some
% factor. These two problems inter infrerence each other. 
% This function ensure that the output values are valid, so the return_flag
% will not be nonzero.

low = 0.03;
high = 0.06; % it is safe to assign high two times bigger than low.

real_min = sqrt(realmin);
real_max = sqrt(realmax);

[ppp pp] = size(cov);

if pp == 1
    if cov <= real_min
        cov = real_min;
    elseif cov >= real_max
        cov = real_max;
    end
    inv_cov = 1 / cov;
    cov_mat = cov;
    det_cov = cov;
    return_flag = 0;
    return;    
end

for i = 1:pp
    if cov(i, i) < real_min
        cov(i, i) =  real_min;
    end
end

if rank(cov) ~= pp
    [eig_vec eig_vl] = eig(cov);
    diag_eig = diag(eig_vl);
    for i = 1:pp
        if isreal(diag_eig(i)) == 0 | diag_eig(i) < 0
            diag_eig(i) = 0;
        end
    end
    sum_eig = sum(diag_eig);
    num_zeros = length(find(diag_eig == 0));
    sort_diag_eig = [diag_eig'; 1:pp];
    sort_diag_eig = sortrows(sort_diag_eig')';
    energy_sum = 0;
    for i = num_zeros + 1 : pp
        energy_sum = energy_sum + sort_diag_eig(1, i);
        if energy_sum < low * sum_eig
            continue;
        elseif energy_sum < high * sum_eig
            sort_diag_eig(1, 1 : i) = energy_sum / i;
            break;
        else
            if i > num_zeros + 1
                sort_diag_eig(1, i) = sort_diag_eig(1, i) - (low * sum_eig - sum(sort_diag_eig(1, num_zeros + 1 : i - 1)));
                sort_diag_eig(1, 1 : i  - 1) = low * sum_eig / (i - 1);
            else
                sort_diag_eig(1, i) = sort_diag_eig(1, i) - low * sum_eig;
                sort_diag_eig(1, 1 : i  - 1) = low * sum_eig / (i - 1);            
            end
            break;
        end    
    end
    diag_eig(sort_diag_eig(2, :)') = sort_diag_eig(1, :)';
    eig_vl = diag(diag_eig);
    inv_cov = real(eig_vec) * (inv(eig_vl)) * real(eig_vec');
    cov_mat = real(eig_vec) * eig_vl * real(eig_vec');
    det_cov = det(eig_vl);
else
    inv_cov = inv(cov);
    cov_mat = cov;
    det_cov = det(cov);
end

return_flag = 0;

if det_cov < real_min
    det_cov = real_min;
    %error('The data should be enlarged, otherwise the inverse function is not accurate.')
    disp('The data should be enlarged, otherwise the inverse function is not accurate');
    %return_flag = 'The data should be enlarged, otherwise the inverse function is not accurate';
elseif det_cov > real_max
    det_cov = real_max;
    disp('The data should be shrinked, otherwise the inverse function is not accurate');
    %error('The data should be shrinked, otherwise the inverse function is not accurate.')
    %return_flag = 'The data should be shrinked, otherwise the inverse function is not accurate';
end



