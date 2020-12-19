%main_simulation.m: test the estimation of simulation for a
%three-compartment model.
%
%
% ========================================================================
%setup by Yuhai Xie 2018
%email: xieyuhai1995@sjtu.edu.cn.
%========================================================================



close all;
clear all;
clc;
load mixedsimdata3SNR10

[A_est,S_est]=CAM(X_mask,3,1,1);

max_a=max(A_est);
[j,g]=sort(max_a);
    
asort(:,1)=A_est(:,g(2));
asort(:,2)=A_est(:,g(1));
asort(:,3)=A_est(:,g(3));

[w1,h1,norm]=NMF(X_mask,3,2000,A_est,S_est);
max_w=max(w1);
[j,g]=sort(max_w);
    
wsort(:,1)=w1(:,g(2));
wsort(:,2)=w1(:,g(1));
wsort(:,3)=w1(:,g(3));
hsort(1,:)=h1(g(2),:);
hsort(2,:)=h1(g(1),:);
hsort(3,:)=h1(g(3),:);

%draw the TC maps
figure;
x=1:size(asort,1);
hl1 = line(x,gtTCp,'LineStyle','--','Color','r');
hold on;
hl2 = line(x,wsort(:,3),'LineStyle','--','Color','g');
hold on;
hl3 = line(x,asort(:,3),'LineStyle','--','Color','b');legend('Ground Truth tissue1','Proposed Method','CAM-CM')

figure;
hl4 = line(x,gtTCf,'LineStyle','--','Color','r');
hold on;
hl5 = line(x,wsort(:,1),'LineStyle','--','Color','g');
hold on;
hl6 = line(x,asort(:,1),'LineStyle','--','Color','b');legend('Ground Truth tissue2','Proposed Method','CAM-CM')

figure;
hl7 = line(x,gtTCs,'LineStyle','--','Color','r');
hold on;
hl8 = line(x,wsort(:,2),'LineStyle','--','Color','g');
hold on;
hl9 = line(x,asort(:,2),'LineStyle','--','Color','b');legend('Ground Truth tissue3','Proposed Method','CAM-CM')

kis=0.1*rand(1);
kos=kis+rand(1);
kif=0.1*rand(1);
kof=kif+rand(1);
initk=[kis,kos;kif,kof];
del_t=0.5;

[eKtrans,eKep,pixelwise_Ktrans]=CM(X_mask,wsort,initk,del_t);
[eKtrans1,eKep1,pixelwise_Ktrans1]=CM(X_mask,asort,initk,del_t);


R1 = double(imread('ROI_fast1.jpg'));  R1=R1(:,:,1)/255; [m n]=size(R1); R1 = ones(m,n)-R1; 
R1(find(R1>0.5))=1; R1(find(R1<=0.5))=0;


R2 = double(imread('ROI_input.jpg'));  R2=R2(:,:,2)/255; R2 = ones(m,n)-R2;
R2(find(R2>0.5))=1; R2(find(R2<=0.5))=0; 

R3 = double(imread('ROI_slow.jpg')); R3=R3(:,:,1)/255; R3 = ones(m,n)-R3;
R3(find(R3>0.5))=1; R3(find(R3<=0.5))=0; 

d=R1+R2+R3;

[r,c]=find(d);

df=zeros(m);
ds=zeros(m);
dp=zeros(m);
df1=zeros(m);
ds1=zeros(m);
dp1=zeros(m);
dhf=zeros(m);
dhs=zeros(m);
dhp=zeros(m);


for i=1:size(r)
    df(r(i),c(i))=pixelwise_Ktrans(1,i);
    ds(r(i),c(i))=pixelwise_Ktrans(2,i);
    dp(r(i),c(i))=pixelwise_Ktrans(3,i);
    df1(r(i),c(i))=pixelwise_Ktrans1(1,i);
    ds1(r(i),c(i))=pixelwise_Ktrans1(2,i);
    dp1(r(i),c(i))=pixelwise_Ktrans1(3,i);
    dhf(r(i),c(i))=hsort(1,i)*eKtrans(1);
    dhs(r(i),c(i))=hsort(2,i)*eKtrans(2);
    dhp(r(i),c(i))=hsort(3,i);

end

%draw the local Ktrans maps
figure;
subplot(2,2,1);imshow(df,[]);colorbar;
subplot(2,2,2);imshow(ds,[]);colorbar;
subplot(2,2,3);imshow(dp,[]);colorbar;
figure;
subplot(2,2,1);imshow(df1,[]);colorbar;
subplot(2,2,2);imshow(ds1,[]);colorbar;
subplot(2,2,3);imshow(dp1,[]);colorbar;
figure;
subplot(2,2,1);imshow(dhf,[]);colorbar;
subplot(2,2,2);imshow(dhs,[]);colorbar;
subplot(2,2,3);imshow(dhp,[]);colorbar;