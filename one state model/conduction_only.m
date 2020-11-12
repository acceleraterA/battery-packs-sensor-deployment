clearvars
close all
clc
%remember to change folder at the end 665/res11 665/res12
%% Battery parameters
Rc = 1.833;%KW-1
%Ru=4.19;%KW-1 
R_e=0.01;%ohm,approximate by equation: Re=V/I^2
Cc = 67;%J/K
Cs=4.5;%J/K
Rcc=2.1;%Al
%% input 
load input_project t I%
% load I_FUDS_newscott_forSPM_2000s_1C25dot67_max3C I 
I=2*I;
I=I';
dt = 0.1;  % sampling period
Ns = 90000;
T_f=[20:10:40];%20,30,40c
%% generate AB matrix
m = input('Please enter the number of the row of the battery system:'); 
n = input('Please enter the number of the column of the battery system:');
if n==1%1*n and m*1 are the same for string of batteries
    n=m;
    m=1;
end
%ru calculator
Dia=0.026; %diameter,m
v=1.515;
% pitch=[1.5 1.25 1.1];
pic=1.5;
% Rcc=Rccc(g);%changing Rcc
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
A = zeros(2*m*n,2*m*n);
B = zeros(2*m*n,2);
x=ones(m*n,Ns).*T_f(1);
%x=[Tc11,Ts11, Tc12,Ts12,...Tc1n,Ts1n, T21,T22....,T2n,...Tm1,Tm2,...,Tmn]'    
for i=2:2:2*m*n
    %%Ts
    %cell to cell conduction
    if m>1
        A(i,i)=-(2/(Ru*Cs)+4/(Rcc*Cs));%all cells except the cells on the boundary have 4 Rccs -4.6546
        A(2,2)=-(2/(Ru*Cs)+2/(Rcc*Cs));
        A(2*n,2*n)=-(2/(Ru*Cs)+2/(Rcc*Cs));
        A(2*(m-1)*n+2,2*(m-1)*n+2)=-(2/(Ru*Cs)+2/(Rcc*Cs));
        A(2*m*n,2*m*n)=-(2/(Ru*Cs)+2/(Rcc*Cs));%4 cells on the corner with 2 Rccs    -2.4323
        if n>2
            for j =4:2:2*(n-1)
                A(j,j)=-(2/(Ru*Cs)+3/(Rcc*Cs));%top row except the corner -3.5435
            end
            for j=2*(m-1)*n+4:2:2*m*n-2
                A(j,j)=-(2/(Ru*Cs)+3/(Rcc*Cs));%bottom row except the corner
            end
        end
        if m>2
            for j=2*n+2:2*n:2*((m-2)*n+1)
                A(j,j)=-(2/(Ru*Cs)+3/(Rcc*Cs));%first column except the corner
            end
            for j=(2*2*n):2*n:(m-1)*2*n
                A(j,j)=-(2/(Ru*Cs)+3/(Rcc*Cs));%last column except corner, cells on the boundary with 3 Rccs
            end
        end
        if mod(i/2,n)~=0
            A(i,i+2)=1/(Rcc*Cs);%right connect rcc except the last column 1.1111
        end
        if mod(i/2,n)~=1
            A(i,i-2)=1/(Rcc*Cs);%left connect rcc except the first column
        end
        if i>2*n+1
            A(i,i-2*n)=1/(Rcc*Cs);%above connect rcc except the first row
        end
        if i<=2*(m-1)*n
            A(i,i+2*n)=1/(Rcc*Cs);%below connect rcc except the last row
        end
        B(i,2)=2/(Ru*Cs);%

    end
end   
A(1:2:2*m*n,:)=[];
A(:,1:2:2*m*n)=[];
B=B(2:2:2*m*n,2);

D=0;

%% enumerate method
k=input('please input the number of sensor:');
for l=4:k     
C1=nchoosek(1:m*n,l);
C_eig=zeros(l,m*n);
N=size(C1,1);
gram_eig_c=zeros(N,1);
gram_tr_c=zeros(N,1);
gram_eig_d=zeros(N,1);
gram_tr_d=zeros(N,1);
% res1={3,1};

for i = 1: N%row of matrix
    C_eig=zeros(l,m*n);
    for j=1:l
    C_eig(j,C1(i,j))=1;
    end
    
    %still discrete
    system_c=ss(A,B,C_eig,D); 
    system_d=c2d(system_c,dt);
    gramians_d = gram(system_d,'o');
    %continuous
       gramians_c = lyap(A',C_eig'*C_eig);
    gram_eig_c(i)= min(eig(gramians_c));%smallest eigenvalue for each combination
    gram_tr_c(i)=trace(gramians_c);
    gram_eig_d(i)= min(eig(gramians_d));%smallest eigenvalue for each combination
    gram_tr_d(i)=trace(gramians_d);
end
res={gram_tr_c;gram_eig_c; gram_tr_d;gram_eig_d};
res1{l,1}=res;
end
save ('number410cond.mat')
%% deployment eig
C1=nchoosek(1:m*n,k);
[m_tr_c,vs_tr_c]=max(res1{k,1}{2,1});
[m_tr_d,vs_tr_d]=max(res1{k,1}{4,1});
C_max=C1(vs_tr_c,:);
c1_c=zeros(k,m*n);
c1_d=zeros(k,m*n);
 l=k;
for i=1:l
        c1_c(i,C1(vs_tr_c,i))=1;
        c1_d(i,C1(vs_tr_d,i))=1;
end%obtain the optimal C 

    %continous
    gramians_c = lyap(A',c1_c'*c1_c);
    eig_c=eig(gramians_c);
    %discrete
    system_c=ss(A,B,c1_d,D); 
    system_d=c2d(system_c,dt);
    gramians_d = gram(system_d,'o');
    eig_d=eig(gramians_d);
    
    %deployment for continous
[x_c,y_c]=num2node(C_max,m,n);
for j=1:length(y_c)
    if y_c(j)==0
        y_c(j)=n;
    end
end
eig123_c(:,:)=zeros(m,n);
for i=1:length(x_c)
    eig123_c(x_c(i),y_c(i))=m_tr_c;
end
eig123_c
%deployment for discrete
% [x_d,y_d]=num2node(C1(vs_tr_d,:),m,n);
% for j=1:length(y_d)
%     if y_d(j)==0
%         y_d(j)=n;
%     end
% end
% eig123_d(:,:)=zeros(m,n);
% for i=1:length(x_d)
%     eig123_d(x_d(i),y_d(i))=m_tr_d;
% end
% eig123_d
%% simulation
%C=zeros(1,m*n);
sys_c=ss(A,B,c1_c,D); 
sys_d = c2d(sys_c,dt);

B_d=sys_d.B; 
A_d=sys_d.A; 
t_end=Ns/1; 
if m>1
    tf=ones((m+1)*n,Ns)*T_f(1);%matrix of tf else
else
    tf=ones(n,Ns)*T_f(1);
end
for w=2:t_end%time step claculation
    x(:,w)=A_d*x(:,w-1)+B_d*I(1,w-1)^2;%Ts and Tc if m>1
end
%% deployment trace

C1=nchoosek(1:m*n,k);
[m_tr,vs_tr]=max(res1{k,1}{1,1});
c_eig=zeros(k,m*n);
for i=1:k
        c_eig(i,C1(vs_tr,i))=1;
end
%observability determined
%     gramians_opt_eig = lyap(A,C_eig'*C_eig);
%     eig1=eig(gramians_opt_eig);
%     if min(eig1)>0
%     disp(['observable'])
%     else
%     disp(['unobservable'])
%     end
%     IsDiagDom(gramians_opt_eig)
%     IsDiagDom(gramians_2)
%     Q=eye(m*n,m*n);
%     gram3=gram(A,Q);
%     IsDiagDom(gram3)
    
    
[x2,y2]=num2node(C1(vs_tr,:),m,n);
for j=1:length(y2)
    if y2(j)==0
        y2(j)=n;
    end
end
trace123(:,:)=zeros(m,n);
for i=1:length(x2)
    trace123(x2(i),y2(i))=m_tr;
end
trace123
% % save ('557cond.mat')
