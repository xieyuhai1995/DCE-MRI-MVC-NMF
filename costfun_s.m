function R=costfun_s(x)
load ('TC.mat')
L=size(TC_s,1);
kis=x(1);
kos=x(2);
for t=1:L
    hs(t)=kis*exp(-(t-1)*del_t*kos);
    Hs(t,1)=hs(t);
    for tt=1:L-t,
        Hs(t+tt,tt+1)=Hs(t,1);
    end
end
R=norm(TC_s-Hs*TC_p);