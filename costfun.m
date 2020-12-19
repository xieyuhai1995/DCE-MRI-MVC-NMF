function R=costfun(x)
load ('TC.mat')
L=size(TC,1);
ki=x(1);
ko=x(2);
for t=1:L
    h(t)=ki*exp(-(t-1)*del_t*ko);
    H(t,1)=h(t);
    for tt=1:L-t,
        H(t+tt,tt+1)=H(t,1);
    end
end
R=norm(TC-H*Cp);