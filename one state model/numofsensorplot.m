clc
clearvars

% ru1=load('445/rcc2.1/v1.5152.11.5.mat');
% % res2=[ru1.res1{1,1};ru1.res1{2,1}{7:8,1}]
% save 445/rcc2.1/v1.5152.11.548.mat res2
ru2=load('number410onestate.mat');
ru3=load('number410cond.mat');
ru4=load('number410conv.mat');
f1=figure()
subplot(2,1,1)
for i=4:1:10
y1(i)=max(ru2.res1{i,1}{1,1});
y2(i)=max(ru2.res1{i,1}{2,1});
y3(i)=max(ru3.res1{i,1}{1,1});
y4(i)=max(ru3.res1{i,1}{2,1});
y5(i)=max(ru4.res1{i,1}{1,1});
y6(i)=max(ru4.res1{i,1}{2,1});

end
plot(4:1:10,y1(4:10),'o-',4:1:10,y3(4:10),'o-',4:1:10,y5(4:10),'o-')
title('(a) tr(W_{o}) versus number of sensors');ylabel('value of tr(W_{o})');
legend('Case 1','Case 2','Case 3')
subplot(2,1,2)
plot(4:1:10,y2(4:10),'o-',4:1:10,y4(4:10),'o-',4:1:10,y6(4:10),'o-')

title('(b) \lambda_{min}(W_{o}) versus number of sensors');xlabel('number of sensors');ylabel('value of \lambda_{min}(W_{o})');
legend('Case 1','Case 2','Case 3')
% xaxis
