%NMF.m: minimun volume constraint non-negative factorization algorithm.
%
%inputs:
%       V: m*n data matrix
%       R: the number of sources 
%       K: the number of iterations
%       w_init: m*K initial combination of vectors of sources
%       h_init: K*n initial distribution of vectors of sources
% ========================================================================
%setup by Yuhai Xie 2018
%email: xieyuhai1995@sjtu.edu.cn.
%========================================================================




function [w,h,dnorm] = NMF(V,R,K,w_init,h_init)

    [m,n]=size(V); 
    tau=0.01;
    
    theta=0.1;
    rho=0.1;
    
    vet_one=ones(1,R);
    vet_delta=ones(1,R)*10;
    vet_delta1=ones(1,n)*10;
    mat_zero=zeros(R-1,R);
    mat_ide=eye(R-1,R-1);
    vet_zero=zeros(1,R-1);
    
    v_aug=[V;vet_delta1];

    mat_c=[vet_one;mat_zero];
    mat_b=[vet_zero;mat_ide];
        
    mat_miu=sum(V,2)/n;
    [~,score]=PCA(V,R-1);
    mat_u=V*score;
       
    w0=w_init;
    h0 =h_init; 
    
    for iter=1:K
      
        alpha=0.1;
        beta=0.1;
        mat_z=mat_c+mat_b*mat_u'*(w0-mat_miu*vet_one);
        inv_z=inv(mat_z);
        w_gradient=(w0*h0-V)*h0'+tau*det(mat_z)^2*mat_u*mat_b'*inv_z';
       
        for i=1:100
            wn=max(0,w0-alpha*w_gradient);d=wn-w0;
            suff_decr=0.5*(sum(sum((V-wn*h0).^2))-sum(sum((V-w0*h0).^2)))+tau*0.5*det(mat_c+mat_b*mat_u'*(wn-mat_miu*vet_one))^2-tau*0.5*det(mat_c+mat_b*mat_u'*(w0-mat_miu*vet_one))^2-theta*alpha*sum(sum(w_gradient.*d))>=0;
            if suff_decr
                alpha=alpha*rho;
            else
                break
            end
        end
        w=max(0,w0-alpha*w_gradient);
        
        w_aug=[w;vet_delta];

        
        h_gradient=w_aug'*(w_aug*h0-v_aug);
        
        for i=1:100
            hn=max(0,h0-beta*h_gradient);d=hn-h0;
            suff_decr=0.5*(sum(sum((v_aug-w_aug*hn).^2))-sum(sum((v_aug-w_aug*h0).^2)))-theta*beta*sum(sum(h_gradient.*d))>=0;
            if suff_decr
                beta=beta*rho;
            else
                break
            end
        end
        h=max(0,h0-beta*h_gradient);
        
        
        dnorm=0.5*sum(sum((V-w*h).^2));
        h0=h;
        w0=w;
    end   
end  


