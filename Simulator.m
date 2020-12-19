clear;clc
close all


% DCE-MRI data simulation, generate the data "simdata.mat"
% make sure you set the current directory to the one containing "Simulator.m"

% % % % % % % % ==============parameter setting ============================================
NumTframe=18;       % number of time frames for Dynamic imaging data
del_t=0.5;          % the sampling interval between time frames (/min)
SNR=10;             % the signal to noise ratio for the simulated data
Ktrans_s = 0.05;    % volume trasfer rate for slow flow
Kep_s =0.5;         % flux rate for slow flow
Ktrans_f = 0.12;    % volume trasfer rate for fast flow
Kep_f = 1.2;        % flux rate for fast flow
% Note: The spatial disbtribution of differnt tissue types are in the figures: ROI_fast.jpg, ROI_slow.jpg, ROI_input.jpg 
% % % % % % % % ==========================================================


%============distribution maps for different tissue types=============
%----------fast volume transfer
R1 = double(imread('ROI_fast1.jpg'));  R1=R1(:,:,1)/255; [m n]=size(R1); R1 = ones(m,n)-R1; 
R1(find(R1>0.5))=1; R1(find(R1<=0.5))=0;
R1(R1==1) = imnoise(R1(R1==1),'gaussian',0,0.05);   %add gaussian noise on the spatial distribution
%----------plasma volume
R2 = double(imread('ROI_input.jpg'));  R2=R2(:,:,2)/255; R2 = ones(m,n)-R2;
R2(find(R2>0.5))=1; R2(find(R2<=0.5))=0; 
R2(R2==1) = imnoise(R2(R2==1),'gaussian',0,0.05);
%--------slow volume transfer
R3 = double(imread('ROI_slow.jpg')); R3=R3(:,:,1)/255; R3 = ones(m,n)-R3;
R3(find(R3>0.5))=1; R3(find(R3<=0.5))=0; 
R3(R3==1) = imnoise(R3(R3==1),'gaussian',0,0.05);

%===========Time activity curves=============
%-------Plasma input function---------
time=0:0.5:15; % mintues
Cp=0.5*(3.99*exp(-1*time)+0.78*exp(-0.0111*time));
%-------Associated fast volume transfer compartment--------- 
Cf = Ktrans_f*conv(Cp,exp(-Kep_f*time));
Cf = Cf(1:length(time));
%-------Associated slow volume transfer compartment--------- 
Cs = Ktrans_s*conv(Cp,exp(-Kep_s*time));
Cs = Cs(1:length(time));

Cp=[0,Cp(1:NumTframe-1)];
Cs=[0,Cs(1:NumTframe-1)];
Cf=[0,Cf(1:NumTframe-1)];


% %-------show Cp, Cs, Cf-----------------
% timeslot=[0:NumTframe-1]*del_t;
% figure; plot(timeslot, Cp, '-b','LineWidth', 4);
% hold on; plot(timeslot, Cs,':g', 'LineWidth', 4);
% hold on; plot(timeslot, Cf,'--r','LineWidth', 4);
% legend ('plasma input','slow flow','fast flow', 'fontsize', 18, 1);
% xlabel('time (min)','fontsize', 18);
% ylabel('Tracer Concentration (mM/L)','fontsize', 18);
% set (gcf,'Position',[300,300,400,300], 'color','w');
% %-------show normalized Cp, Cs, Cf-----------------
% figure; plot(timeslot, Cp/sum(Cp), '-b','LineWidth', 4);
% hold on; plot(timeslot, Cs/sum(Cs),':g', 'LineWidth', 4);
% hold on; plot(timeslot, Cf/sum(Cf),'--r','LineWidth', 4);
% legend ('plasma input','slow flow','fast flow', 'fontsize', 18, 1);
% xlabel('time (min)','fontsize', 18);
% ylabel('Normalized TC','fontsize', 18);
% set(gca,'fontsize',18)


%========Data generation=========
sumR=R1(:)+R2(:)+R3(:)+eps;
S = [R1(:)'./sumR'; R2(:)'./sumR'; R3(:)'./sumR'];
A = [Cs' Cp' Cf']; 
X = A*S; 
%------------loading the mask------------------
index=find(sum(S,1)>0);
index=index';
N=length(index); % N is the number of the pixels in masked ROI
M=NumTframe;

X_mask=X(:,index);
Rn=eye(M);
alpha=0.9;
for i=1:M
    for j=1:M
        Rn(i,j)=alpha^(abs(i-j));
    end
end
R=chol(Rn);

for i=1:N,
    sigma=std(X_mask(:,i))*10^(-SNR/10);
    X_mask(:,i)=X_mask(:,i)+sigma*R'*randn(M,1);
end
X_mask(find(X_mask<0))=eps;
origX=X(:,index);
gtS=S(:,index);  gtSs=gtS(1,:); gtSp=gtS(2,:); gtSf=gtS(3,:);
gtTCp=Cp';gtTCf=Cf';gtTCs=Cs';
save mixedsimdata3SNR10 gtTCf gtTCs gtTCp gtSs gtSp gtSf del_t Ktrans_f Kep_f Ktrans_s Kep_s X_mask 