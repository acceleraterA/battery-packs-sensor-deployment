
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
%% input 
load input_project I t %I = 5A; t=1:90000
% load I_FUDS_newscott_forSPM_2000s_1C25dot67_max3C I 
I=2*I;
dt = 0.1;  % sampling period
Ns = length(t);%90000
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
pic=1.5;
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
        A(i,i)=-2/(Ru*Cs);%all cells except the cells on the boundary have 4 Rccs -4.6546
        A(2,2)=-2/(Ru*Cs);
        A(2*n,2*n)=-2/(Ru*Cs);
        A(2*(m-1)*n+2,2*(m-1)*n+2)=-2/(Ru*Cs);
        A(2*m*n,2*m*n)=-2/(Ru*Cs);%4 cells on the corner with 2 Rccs    -2.4323
        if n>2
            for j =4:2:2*(n-1)
                A(j,j)=-2/(Ru*Cs);%top row except the corner -3.5435
            end
            for j=2*(m-1)*n+4:2:2*m*n-2
                A(j,j)=-2/(Ru*Cs);%bottom row except the corner
            end
        end
        if m>2
            for j=2*n+2:2*n:2*((m-2)*n+1)
                A(j,j)=-2/(Ru*Cs);%first column except the corner
            end
            for j=(2*2*n):2*n:(m-1)*2*n
                A(j,j)=-2/(Ru*Cs);%last column except corner, cells on the boundary with 3 Rccs
            end
        end
%         if mod(i/2,n)~=0
%             A(i,i+2)=1/(Rcc*Cs);%right connect rcc except the last column 1.1111
%         end
%         if mod(i/2,n)~=1
%             A(i,i-2)=1/(Rcc*Cs);%left connect rcc except the first column
%         end
%         if i>2*n+1
%             A(i,i-2*n)=1/(Rcc*Cs);%above connect rcc except the first row
%         end
%         if i<=2*(m-1)*n
%             A(i,i+2*n)=1/(Rcc*Cs);%below connect rcc except the last row
%         end
        %coolant and battery surface convection
        if i>2*n &&mod(i/2,n)~=1 %start at second row cell and except for first column which tf = tin
            for j =2:2:2*(i/2-(ceil(i/2/n)-1)*n-1)%ceil(i/2/n)=#column
                if i>2*n && i< 2*(m-1)*n +1&& m>2 %except the first and last row
                    A(i,i-j)=A(i,i-j)+2*(1-2/(Ru*Cf))^(j/2-1)*(1/(Ru^2*Cf*Cs));
                    A(i,i-j-2*n)=A(i,i-j-2*n)+(1-2/(Ru*Cf))^(j/2-1)*(1/(Ru^2*Cf*Cs));
                    A(i,i-j+2*n)=A(i,i-j+2*n)+(1-2/(Ru*Cf))^(j/2-1)*(1/(Ru^2*Cf*Cs));
                elseif i> 2*(m-1)*n+2%last row
                    A(i,i-j)=A(i,i-j)+(1-2/(Ru*Cf))^(j/2-1)*(1/(Ru^2*Cf*Cs))+(1-1/(Ru*Cf))^(j/2-1)*(1/(Ru^2*Cf*Cs));%current row
                    A(i,i-j-2*n)=A(i,i-j-2*n)+(1-2/(Ru*Cf))^(j/2-1)*(1/(Ru^2*Cf*Cs));%row above
                end
            end
        elseif i<=2*n %first row except first column
            for j=2:2:i-2
                A(i,j)=A(i,j)+(1-2/(Ru*Cf))^(j/2-1)*(1/(Ru^2*Cf*Cs))+(1-1/(Ru*Cf))^(j/2-1)*(1/(Ru^2*Cf*Cs));%current row
                A(i,j+2*n)=A(i,j+2*n)+(1-2/(Ru*Cf))^(j/2-1)*(1/(Ru^2*Cf*Cs));%row below
            end
        end
        B(i,1)=R_e/Cs;
        B(i,2)=2/(Ru*Cs)*(1-2/(Ru*Cf))^(i/2-(ceil(i/2/n)-1)*n-1);         

    else% for a string of battries
        %conduction between each battries
        A(i,i)=-(1/(Ru*Cs));
        A(2,2)=-(1/(Ru*Cs));
        A(2*n,2*n)=-(1/(Ru*Cs));%first and last cells have one Rcc
        
        %coolant and battery surface convection
        for j = 2:2:i-2
            A(i,j)=A(i,j)+(1-1/(Ru*Cf))^((i-j)/2-1)/(Ru^2*Cf*Cs);
        end
        B(i,1)=R_e/Cs;
        B(i,2)=1/(Ru*Cs)*(1-1/(Ru*Cf))^(i/2-1);

    end
end
A(1:2:2*m*n,:)=[];
A(:,1:2:2*m*n)=[];
B=B(2:2:2*m*n,1);
D=0;
k=input('please input the number of sensor:');
% k=4;


%% enumerate method

for l=4:k     
C1=nchoosek(1:m*n,l);
C_eig=zeros(l,m*n);
N=size(C1,1);
gram_eig_c=zeros(N,1);
gram_tr_c=zeros(N,1);
% gram_eig_d=zeros(N,1);
% gram_tr_d=zeros(N,1);

for i = 1: N%row of matrix
    C_eig=zeros(l,m*n);
    for j=1:l
    C_eig(j,C1(i,j))=1;
    end
%     system_c=ss(A,B,C_eig,D); 
%     system_d=c2d(system_c,dt);
% 
%     gramians_c = gram(system_c,'o');
    %continuous
    gramians_c = lyap(A',C_eig'*C_eig);   
    gram_eig_c(i)= min(eig(gramians_c));%smallest eigenvalue for each combination
    gram_tr_c(i)=trace(gramians_c);
%     gram_eig_d(i)= min(eig(gramians_d));%smallest eigenvalue for each combination
%     gram_tr_d(i)=trace(gramians_d);
end
res={gram_tr_c;gram_eig_c};
res1{l,1}=res;
end
% mkdir results
% save ('results/conv.mat')
save ('number410conv.mat')
%% deploy eig
C1=nchoosek(1:m*n,k);
[m_tr_c,vs_tr_c]=max(res1{k,1}{2,1});
c1_c=zeros(k,m*n);
l=k;
for i=1:l
        c1_c(i,C1(vs_tr_c,i))=1;
end%obtain the optimal C 
    %continous
%         system_c=ss(A,B,c1_c,D); 
%     system_d=c2d(system_c,dt);
% 
%     gramians_c = gram(system_c,'o');
    gramians_c = lyap(A',c1_c'*c1_c);
    eig_c=eig(gramians_c);
    %deployment for continous
[x_c,y_c]=num2node(C1(vs_tr_c,:),m,n);
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
%% simulation
C=zeros(1,m*n);
sys_c=ss(A,B,C,D); 
sys_d = c2d(sys_c,dt);

B_d=sys_d.B; 
A_d=sys_d.A; 
t_end=Ns/1; 
for w=2:t_end%time step claculation
    x(:,w)=A_d(:,:)*x(:,w-1)+B_d(:,:)*I(w-1)^2;%Ts and Tc if m>1
end

%% deployment trace
C1=nchoosek(1:m*n,k);
[m_tr,vs_tr]=max(res1{k,1}{1,1});
c_eig=zeros(k,m*n);
for i=1:k
        c_eig(i,C1(vs_tr,i))=1;
end

[x2,y2]=num2node(C1(vs_tr,:),m,n);
for j=1:length(y2)
    if y2(j)==0
        y2(j)=n;
    end
end
eig4(:,:)=zeros(m,n);
for i=1:length(x2)
    eig4(x2(i),y2(i))=m_tr;
end
eig4
% % save ('557conv.mat')