function [normindic, indic, estmu, estcov, estpp, id_record, return_flag] = SL_EM(y, estmu, estcov, estpp, min_num_point, flag_allow_vanish)

%%%%%%%% Attention: This function can not be used with data dimension
%%%%%%%% bigger than 1 and MDL criterion
% y is the data, each column is a sample
% estmu is the cluster center, each column is a center
% estcov is the covariance matrix, 
% estpp is the mix proportion of each cluster
% min_num_point is the minimum number of samples in a cluster
% flag_allow_vanish, 1: allow cluster to vanish in the EM process, 0:
% doesn't allow cluster to vanish, so return err
% id_record contains the id of clusters that remains after EM.
% return_flag: 0, everything is ok, 1, cluster vanished and return when
% cluster vanish is not allowed, 2, not converged return.

real_min = sqrt(realmin);
real_max = sqrt(realmax);

[dimens, npoints] = size(y);
[dimens, k] = size(estmu);
id_record = 1:k;

indic = zeros(k, npoints);
semi_indic = zeros(k, npoints);
flag = 1;
while flag
    for i=1:k
        temp_c = zeros(1, dimens);
        temp_c(find(diag(estcov(:, :, i))' < real_min)) = real_min;
        estcov(:, :, i) = estcov(:, :, i) + diag(temp_c);        
        semi_indic(i,:) = multinorm(y, estmu(:,i), estcov(:,:,i));
        indic(i,:) = semi_indic(i, :) * estpp(i);
    end
    if length(find(sum(indic, 1) < real_min)) > 0
        estcov = 2 * estcov;
    else
        flag = 0;
    end
end
    
cont=1;        % auxiliary variable of the inner loop
countstep = 0;
th = 0.001;
pre_normindic = ones(k, npoints);
normindic = ones(k, npoints);
err_team = [];
ML_team = [];
while(cont)
    % E step
    clear semi_indic indic normindic;
    for i=1:k
        semi_indic(i,:) = multinorm(y, estmu(:,i), estcov(:,:,i));
        indic(i,:) = semi_indic(i,:) * estpp(i);
    end
    t = find(sum(indic, 1) < real_min);
    num_t = length(t);
    if num_t > 0
        for i = 1:num_t
            if min(indic(:, t(i))) == max(indic(:, t(i)))
                indic(:, t(i)) = real_min * ones(k, 1) / k;
            else
                indic(:, t(i)) = indic(:, t(i)) / sum(indic(:, t(i)), 1) * real_min ;
            end
        end
    end
    normindic = indic ./ kron(ones(k,1), sum(indic, 1));
    % M step
    estpp = sum(normindic') / npoints;
    id_vanish = find(estpp * npoints < min_num_point);
    if length(id_vanish) == 0
        clear estmu estcov;
        for i = 1:k
            normalize = 1/sum(normindic(i, :));
            estmu(:, i) = y * normindic(i, :)' * normalize;
            aux = kron(normindic(i, :), ones(dimens, 1)) .* (y - estmu(:, i) * ones(1, npoints));
            estcov(:, :, i) = normalize * (aux * (y - estmu(:, i) * ones(1, npoints))');
            temp_c = zeros(1, dimens);
            temp_c(find(diag(estcov(:, :, i))' < real_min)) = real_min;
            estcov(:, :, i) = estcov(:, :, i) + diag(temp_c);
        end
    else
        if flag_allow_vanish == 0
            return_flag = 1;
            normindic = 0; 
            indic = 0; 
            estmu = 0; 
            estcov = 0;
            estpp = 0;
            id_record = 0;
            return;
        end
        [min_estpp id_vanish] = min(estpp);        
        id_record(id_vanish) = [];
        remain_id = 1:k;
        remain_id(id_vanish) = [];
        estmu = estmu(:, remain_id);
        estcov = estcov(:, :, remain_id);
        estpp = estpp(remain_id);
        estpp = estpp / sum(estpp);
        k = k - 1; 
        %disp([num2str(length(id_vanish)) ' clusters vanished.' num2str(k) ' clusters remain.']);
    end    
    
    [num_com num_d] = size(pre_normindic);
    countstep = countstep + 1;  
    if num_com == k
        err = mean(sum(abs(normindic - pre_normindic), 1)); 
        err_team = [err_team err];
        ML = sum(log(sum(indic, 1)));
        ML_team = [ML_team ML];
        if err < th
           cont = 0;
           err;
           ML_team(max(1, length(ML_team) - 2) : length(ML_team));
           estpp;
           countstep;           
        elseif rem(countstep, 1000) == 0
            err_team(length(err_team) - 10 : length(err_team));
            ML_team(length(ML_team) - 10 : length(ML_team));
            disp('can not converge');
            % s = input('1 for continue, 2 for stop:  ');
            s = 2;
            if s == 2
                return_flag = 2;
                return;
            end
        end
    else
        countstep = 0;
    end
    pre_normindic = normindic;
end % this end is of the inner loop: "while(cont)"
return_flag = 0;