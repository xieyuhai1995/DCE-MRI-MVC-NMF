function [eA,cornerind]=measure_conv(X,N)

[M,L]=size(X);
total_index=nchoosek(1:L,N);
comb=size(total_index,1);

error = zeros(comb,1);
for p=1:comb
    A = X(:,total_index(p,:));
    Others = X(:,setdiff(1:L,total_index(p,:)));
    Ae = [1e-5*A; ones(1,N)]; 
    Oe = [1e-5*Others; ones(1,L-N)];
    alpha = zeros(N,L-N);
    for i=1:L-N
        opt = optimset('Display', 'off');
        alpha(:,i)=lsqnonneg(Ae,Oe(:,i),[],opt);
    end
    error(p) = norm(Others - A*alpha,'fro')^2;
end
[val ind]=min(error);

eA=X(:,total_index(ind,:));
cornerind=total_index(ind,:);