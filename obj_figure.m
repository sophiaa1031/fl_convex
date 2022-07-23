%% 目标函数的值以及变化关系
clc;clear;close all

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
% f = (1e9:1e8:2e9);
f = 1e9;
distance = 25;%[2,50]
v = 30;

g=10;
loss_function = zeros(1,g);
T_max = 0.6;
erequirement = 0.005;
p = 0.001;
q = (0.04:0.01:1);  
% p = (0.001:0.001:0.01);
% q = 0.04;
%%
obj = exp((I+bkn0)*inv_pos(p)./(distance.^(-gamma))*(2.^(z./(b* q.*(T_max-itr*D*c./f)))-1))./q;

 %%
figure;
plot(q,obj,'b-o','Linewidth',2,'MarkerSize',8);
grid on
% xlabel('maximum latency');
xlabel('q');
ylabel('loss function upper bound');