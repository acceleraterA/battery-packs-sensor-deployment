clearvars
close all
clc
%remember to change folder at the end 665/res11 665/res12
%% Battery parameters
dote=16;
Rc = 1.833;%KW-1
%Ru=4.19;%KW-1 
R_e=0.01;%ohm,approximate by equation: Re=V/I^2
Cc = 67;%J/K
Cs=4.5;%J/K
% rcc=linspace(0.2,2,10);%Al
% rcc=linspace(0.1,10,dote);
% Rcc=2;
%Cf=6.1; 
% T_f=[20:10:40];%20,30,40c

%% generate matrix
%m = input('Please enter the number of the row of the battery system:'); 
%n = input('Please enter the number of the column of the battery system:');
m=4;
n=4;
%ru calculator
Dia=0.026; %diameter,m

gram_eig_one2=zeros(dote,1);
gram_tr_one2=zeros(dote,1);
gram_eig_cond=zeros(dote,1);
gram_tr_cond=zeros(dote,1);
gram_eig_conv=zeros(dote,1);
gram_tr_conv=zeros(dote,1);

k=6;
ki=[1 3 8 10 12 14];
C_eig=zeros(k,m*n);
for i=1:k
    C_eig(i,ki(i))=1;
end
D=0;
v=1.5;
pic=1.5;
for p=1:dote
rcc=linspace(0.1,5,dote);
Rcc=rcc(p);%changing Rcc
%ru calculator
vis=1.562*1e-5;%Kinematic viscosity(Air),(m^2/s)
C2=[0.7 0.8 0.86 0.9 0.92 0.94 0.95 0.95 0.96 0.97];%number of row 1-10
pr=0.7296;
tk=0.02551;%Thermal conductivity(Air)(w/(m*K))
L=0.065;%length,m
doa=1.184;%density of air(25C)(kg/m^3)
cp=1007;%specific heat of airJ/(kg*K)
ST=pic*Dia;%changing 
vmax=ST/(ST-Dia)*v;
Re1= vmax*Dia/vis;%renolds number
Nu=C2(m)*0.27*(Re1^0.63)*(pr^0.36);

Vdot=ST*L*v/2;
NTU=pi*Nu*tk*L/(doa*Vdot*cp);
Ru=1/(doa*Vdot*cp*(1-exp(-NTU)));
Cf=doa*Vdot*cp;

[A_one,B_one] = A_onestate_make(Rcc,Ru,Cs,Cf,R_e,m,n);
[A_cond,B_cond] = A_cond_make(Rcc,Ru,Cs,m,n);
% [A_conv,B_conv] = A_conv_make(Ru,Cs,Cf,R_e,m,n);

    gram_eig_one2(p,1)= min(eig(lyap(A_one',C_eig'*C_eig)));%smallest eigenvalue for each combination
    gram_tr_one2(p,1)=trace(lyap(A_one',C_eig'*C_eig));
    
    gram_eig_cond(p,1)= min(eig(lyap(A_cond',C_eig'*C_eig)));%smallest eigenvalue for each combination
    gram_tr_cond(p,1)=trace(lyap(A_cond',C_eig'*C_eig));
    
%     gram_eig_conv(p,1)= min(eig(lyap(A_conv',C_eig'*C_eig)));%smallest eigenvalue for each combination
%     gram_tr_conv(p,1)=trace(lyap(A_conv',C_eig'*C_eig));  
%             
end
p=1:16;
load st_case1_conv.mat
figure()
subplot(2,1,1)
plot(rcc(p),gram_tr_one2(p,1),rcc(p),gram_tr_cond(p,1))
xlabel('R_{cc}');ylabel('tr(W_{o})');
title('(a)');%with sensors located at [1 3 8 10 12 14], (S_{T}/D=1.5
legend('case1','case2');
subplot(2,1,2)
plot(rcc(p),gram_eig_one2(p,1),rcc(p),gram_eig_cond(p,1))
xlabel('R_{cc}');ylabel('\lambda_{min}(W_o)');
title('(b)');% R_{cc} vs \lambda_{min}(W_{o})with sensors located at [1 3 8 10 12 14], (S_{T}/D=1.5)
legend('case1','case2');
% subplot(2,2,3)
% plot(pitch(p),gram_tr_one(p,1),pitch(p),gram_tr_conv(p,1))
% xlabel('S_{T}/D');ylabel('tr(W_{o})');
% title('S_{T}/D vs Tr(W_{o}) with sensors located at [1 3 4 5 7 8 9 10 12 13 14 16], R_{cc}=2)');
% legend('case1','convection');
% subplot(2,2,4)
% plot(pitch(p),gram_eig_one(p,1),pitch(p),gram_eig_conv(p,1))
% xlabel('S_{T}/D');ylabel('\lambda_{min}(W_o)');
% title('S_{T}/D vs \lambda_{min}(W_{o}) with sensors located at [1 3 4 5 7 8 9 10 12 13 14 16], R_{cc}=2)');
% legend('case1','convection');


% p=1:dote;
% figure()
% plot(pitch(p), CF(p),pitch(p),RU(p))
% legend('C_{f}','R_{u}');
% xlabel('S_{T}');
% title('changing of S_{T} on the value of C_{f} and R_{u}');
