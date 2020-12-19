function [R,eigV]=PCA(A,aimdim)
% A is input, and each row is a data point
% R is output, which is the dimension decreased data, each row is a data
% point
num=size(A,1);  %row number
M=size(A,2);    %M is dimension num
      %aimdim is the dimension you want to reduce to

% PCA must do on zero mean data
for i=1:num,
    B(i,:)=A(i,:)-mean(A,1);
end

C=cov(B);
[V,D]=eig(C);
% D is diagnal matrix of eigenvalues of C, D is corresponding eigenvectors
% A*V=V*D

% if you need dimension reduction, then, just do something on D

newD=D(M-aimdim+1:M,M-aimdim+1:M);
newV=V(:,M-aimdim+1:M);
R=B*newV;   %output dimension-descreased data
eigV=newV;
