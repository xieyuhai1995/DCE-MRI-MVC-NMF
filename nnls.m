function S=nnls(X,A)

[M L]=size(X);  %M is dimension, L is number of points
N=size(A,2);
S = zeros(N,L);
for i=1:L
    opt = optimset('Display', 'off');
    S(:,i)=lsqnonneg(A,X(:,i),[],opt);
end
