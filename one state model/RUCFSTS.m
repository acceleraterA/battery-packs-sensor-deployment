Dia=0.026; %diameter,m
%Changing velocity
% veloc=[1 1.515 2 3];
v=1.515;
m=4;

% pitch=[1.5 1.25 1.1];
% k=input("Please enter the maximum sensor number:");%enumerate method
% for g=1:length(pitch)
%     pic=1.1;%changing velocity
dote=20;
pic=linspace(1.1,10,dote);

% Rcc=Rccc(g);%changing Rcc
%ru calculator
vis=1.846*1e-5;%Kinematic viscosity(Air),(m^2/s)
C2=[0.7 0.8 0.86 0.9 0.92 0.94 0.95 0.95 0.96 0.97];%number of row 1-10
pr=0.707;
tk=0.0263;%Thermal conductivity(Air)(w/(m*K))
L=0.065;%length,m
doa=1.161;%density of air(25C)(kg/m^3)
cp=1007;%specific heat of air J/(kg*K)
rus=zeros(dote,1);
cf=zeros(dote,1);
Vmax=zeros(dote,1);
% vdot=v*ST*L;
R_e1=zeros(dote,1);
area=pi*Dia*L;

for i=1:dote
ST=pic(i)*Dia;%changing 
vmax=ST/(ST-Dia)*v;
Re1= vmax*Dia/vis;%renolds number
Nu=0.89*0.35*(Re1^0.6)*(pr^0.36);
% Nu=0.27*Re1^0.63*(pr^0.36)*(Pr/Prw)^0.25;
NTU=pi*Nu*tk/(doa*v*cp*ST);

Ru=1/(doa*ST*L*v*cp*(1-exp(-NTU)));
Cf=doa*v*ST*L*cp;
rus(i)=Ru;
cf(i)=Cf;
Vmax(i)=vmax;
R_e1(i)=Re1;
end
q=1:dote;
figure()
plot(pic(q),rus(q))
% title('ru and cf when changing ST');
xlabel('S_{T}/D');
ylabel('R_{u}');
% legend('Ru','C);
% figure()
% plot(st(q),v
save rus
% figure()
% i=1:dote;
% plot(pic(i),Ru(i),pic(i),rus(i)) 
% legend('align','staggered');
% xlabel('ST/D');ylabel('ru');title('comparing Ru with align and staggered');
