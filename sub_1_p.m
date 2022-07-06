%% clc;clear;close

m=50; %Users' total number
iterationNumber=10;    %Maximun iteration times
hk=0.000004;    %Channel gain 
I=3.7*10^(-8);  %Interference
bkn0=10^(-14);  %Noise = B * N0
b= 1e6; %bandwidth
gamma=3.76;
p0=0.01;  %uplink power
itr = 2;%local iterations
Dia = 15; %diameter of the BS
z = 28.1; %model size
T_max = 35*ones(1,m); %maximum latency in one round
erequirement=0.7*ones(1,m); %% energy requirement
kapa = 1e-28;%processor coefficient
D = 128;%sample size
c = 1.68e7;%cpu bits

f=normrnd(10^9,10^8,1,m);  %CPU frequency
distance = normrnd(7,2,1,m); %distance to the BS
v = normrnd(13,5,1,m); %velocity

% f= 7e8+(1e9-7e8)*rand(1,m);
% distance = 1+12.*rand(1,m);
% v = 30.*rand(1,m);

q = rand(1,m);%RB assignment
p = ones(1,m)*0.1; % 决策变量功率p初始化
lambda = ones(m,1)*0.05;   % 拉格朗日乘子lambda初始化 
real_ite=1; % 当前迭代次数
G_v = [];  
lambda_v = [];  
L_v = [];  % 拉格朗日函数
obj_v = []; 
tol = 10^(-4); % 判断收敛的阈值
t_lambda=0.01;  % 拉格朗日乘子的更新步长
t_x=0.01;      % 决策变量 x 的更新步长
max_iteration = 3e4;  % 最大迭代次数

for i = 1 : max_iteration
    %% 迭代更新功率P
    for k=1:m
        H_1 = (I+bkn0)*(2^(z/(b*q(k)*(T_max(k)-itr*D*c/f(k)))-1))/...
            (distance(k)^(-gamma));
        H_2 = z/q(k)/b;
        H_3 = hk/(I+bkn0);
        H_4 = -H_2*H_3*p(k)/log(2)/(1+H_3*p(k))/(log(1+H_3*p(k))/log(2))^2+...
            H_2/(log(1+H_3*p(k))/log(2));
        part_p = -H_1*exp(H_1/p(k))/(q(k)*(p(k))^2)+lambda(k)*H_4;
        
        p(k) = p(k)-t_x*part_p;
    end
    %% 迭代更新拉格朗日乘子 lambda
    for k=1:m
        Grad_f(k) = p(k)*z*inv_pos(q(k)*b*log(1+p(k)*hk/(I+bkn0))/log(2))+itr*kapa*D*c*f(k)^2 - erequirement(k);
    end
    lambda = max( lambda + t_lambda * Grad_f' , zeros( m,1) );

    for k=1:m
    Term1(k) = 2^(z/(b*q(k)*(T_max(k)-itr*D*c/f(k)))-1);
    end
    obj = 0;
    for k=1:m
        obj = obj + exp((I+bkn0)*Term1(k)*inv_pos(p(k))/(distance(k)^(-gamma)))/q(k);
    end
    
    % 拉格朗日函数
    L = obj+sum(Grad_f * lambda);
    % 记录历史更新结果
    G_v = [G_v Grad_f];%乘子的梯度
    lambda_v = [lambda_v lambda];%乘子的值
    obj_v = [obj_v obj]; % obj_v append
    L_v = [L_v L]; % L_v append
    %% 终止条件判断
    if (real_ite>2) && ( abs( obj - L ) < tol )  % k == 5000 %k>3&&(f_v(k)-f_v(k-1)<tol)
        break;
    end
    real_ite=real_ite+1;
end


if i == max_iteration
    disp(['超过最大迭代次数！']);
end
disp(['迭代次数为：',num2str(real_ite),'次']);
disp(['最优值点为：']);
p

%% 查看收敛性
figure(1)
plot( 2:length(L_v),L_v(2:length(L_v)), '-ro', 'linewidth',2)
grid on
xlabel('Iteration round');
ylabel('Lagrange function');

figure(2)
plot( 2:length(L_v), obj_v(2:length(obj_v)) - L_v(2:length(L_v)) , '-b*', 'linewidth',2)
grid on
xlabel('Iteration round');
ylabel('Objective - Lagrange function');

%% p与f，d关系
figure(3)
pfd=[(p*100)',(f/1e8)',distance']
bar(1:10,pfd(1:10,:))
legend("p*10^(-2)","f*10^8","distance")