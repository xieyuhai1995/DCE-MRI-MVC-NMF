function R=costfun_f(x)
load ('TC.mat')
L=size(TC_f,1);
kif=x(1);
kof=x(2);
for t=1:L
    hf(t)=kif*exp(-(t-1)*del_t*kof);
    Hf(t,1)=hf(t);
    for tt=1:L-t,
        Hf(t+tt,tt+1)=Hf(t,1);
    end
end
R=norm(TC_f-Hf*TC_p);