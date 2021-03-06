%% 寻找合适固定参数
clc;
%% 
clear;close all

m=10; %Users' total number
N = m*0.5; %the number of resource blocks
hk=4e-6;    %Channel gain
I=0;  %Interference
bkn0=1e-13;  %Noise = B * N0
b= 1e7; %bandwidth
gamma=2.5;
Dia = 15; %diameter of the BS
z = 1e6; %model size 100Mb
kapa = 1e-29;%processor coefficient
D = z;%sample size
c = 40;%cpu bits
T = 200;%T=k*I 所有训练轮数
itr =5;
f = 1e9;%[1,3]
distance = 25;%[2,50]
T_max = 0.8;%[0.6,0.86]
erequirement = 0.005;%[2.5e-3,5.5]

% v = (0.001:0.001:0.01);
% q = 0.04;
p = 0.005;
v = (0.04:0.01:1);  



%%
k = 1;
for q=0.04:0.01:1
    p_u(k) = z*p/(b* q*log2(1+p*hk*(distance^(-gamma))/(I+bkn0)));
    p_c(k) = itr*kapa*D*c*(f)^2;
    t_u(k) = z/(b* q*log2(1+p*hk*(distance^(-gamma))/(I+bkn0)));
    t_c(k) = itr*D*c/f;
    k = k+1;
end
figure;
plot(v,p_u,'b-o',v,p_c,'r-o','Linewidth',2,'MarkerSize',8);
grid on
xlabel('p');
ylabel('power');
legend('uploading','computing')
figure;
plot(v,t_u,'b-o',v,t_c,'r-o','Linewidth',2,'MarkerSize',8);
grid on
xlabel('p');
ylabel('latency');
legend('uploading','computing')

%%
% p_u = z*p/(b* q*log2(1+p*hk*(distance^(-gamma))/(I+bkn0)));
% p_c = itr*kapa*D*c*(f)^2;
% t_u = z/(b* q*log2(1+p*hk*(distance^(-gamma))/(I+bkn0)));
% t_c = itr*D*c/f;
% display(['uploading power is: ',num2str(p_u),', computing power is: ',num2str(p_c)]);
% display(['uploading time is: ',num2str(t_u),', computing time is: ',num2str(t_c)]);
