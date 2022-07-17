clc;clear;close all

m=50; %Users' total number
N = m*0.5; %the number of resource blocks
iterationNumber=10;    %Maximun iteration times
hk=0.000004;    %Channel gain
I=3.7*10^(-8);  %Interference
bkn0=10^(-14);  %Noise = B * N0
b= 1e6; %bandwidth
gamma=3.76;
Dia = 15; %diameter of the BS
z = 1e4; %model size 10Mb
T_max = 35*ones(1,m); %maximum latency in one round
erequirement=0.7*ones(1,m); %% energy requirement
kapa = 1e-28;%processor coefficient
D = 128;%sample size
c = 1.68e7;%cpu bits
T = 200;%T=k*I 所有训练轮数

% f=normrnd(10^9,10^8,1,m);  %CPU frequency
% distance = normrnd(7,2,1,m); %distance to the BS
% v = normrnd(13,5,1,m); %velocity

f= 7e8+(1e9-7e8)*rand(1,m);
% distance = 1+12.*rand(1,m);
% f= 1e8*ones(1,m);
distance = 1+50.*rand(1,m);
v = 30.*rand(1,m);

% p = ones(1,m)*0.1; % 功率p初始化
p = ones(1,m)*0.5; % 功率p初始化
q = ones(1,m)*0.1;%RB assignment初始化
itr = 1;%local iterations初始化

p_int = p; % 功率p初始化
q_int = q;%RB assignment初始化
itr_int = itr;%local iterations初始化
 
epsilon1 = ones(1,m)*0.0001; %p的收敛条件
epsilon2 = ones(1,m)*0.0001; %q的收敛条件
p_last = zeros(1,m);
q_last = zeros(1,m);
round_total = 0;

 while(sum(abs(p_last-p)>=epsilon1)>0 && sum(abs(q_last-q)>=epsilon2)>0)
%     p_last = p.copy();
    p_last = p;
    sub_1_p; %更新p

    q_last =q;
    sub_2_q;%更新q
    
    sub_3_I;
    
    round_total = round_total+1;
end
disp(['********迭代次数为：***********',num2str(round_total),'次']);
disp(['最优值点为：']);
p
q
itr
