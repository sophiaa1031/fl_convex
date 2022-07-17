% clc;clear;close all

% m=50; %Users' total number
% N = m*0.5; %the number of resource blocks
% iterationNumber=10;    %Maximun iteration times
% hk=0.000004;    %Channel gain
% I=3.7*10^(-8);  %Interference
% bkn0=10^(-14);  %Noise = B * N0
% b= 1e6; %bandwidth
% gamma=3.76;
% Dia = 15; %diameter of the BS
% z = 28.1; %model size
% T_max = 35; %maximum latency in one round
% erequirement=70*ones(1,m); %% energy requirement
% kapa = 1e-28;%processor coefficient
% D = 128;%sample size
% c = 1.68e7;%cpu bits
% T = 200;%T=k*I 所有训练论数
% 
% % f=normrnd(10^9,10^8,1,m);  %CPU frequency
% % distance = normrnd(7,2,1,m); %distance to the BS
% % v = normrnd(13,5,1,m); %velocity
% 
% f= 7e8+(1e9-7e8)*rand(1,m);
% distance = 1+12.*rand(1,m);
% v = 30.*rand(1,m);
% p = rand(1,m); % 功率p初始化
% q = rand(1,m);%RB assignment
% itr = 2;%local iterations
part_itr =zeros(1,m);%I求导暂时存值
lambda = ones(m,1)*0.05;   % 拉格朗日乘子lambda初始化 
real_ite=1; % 当前迭代次数
G_v = [];  
lambda_v = [];  
L_v = [];  % 拉格朗日函数
obj_v = []; 
tol = 10^(-4); % 判断收敛的阈值
t_lambda=0.01;  % 拉格朗日乘子的更新步长
t_x=0.01;      % 决策变量 x 的更新步长
max_iteration = 3e3;  % 最大迭代次数
for i = 1 : max_iteration
    %% 迭代更新本地轮数I
    for k=1:m
        v1 = (I+bkn0)/(p(k)*distance(k)^(-gamma));
        v2 = 2^(z*f((k))/(b*q(k)*D*c));
        v3 = T_max(k) *f(k)/(D*c);
        v4 = p(k)*z*T*inv_pos(q(k)*itr^2*b*log(1+p(k)*hk/(I+bkn0))/log(2));
        part_itr(k) = exp(2^(v2/(v3-itr)-1)*v1)*(4*v3^2-8*v3*itr+itr*(4*itr+2^(v2/(v3-itr))*v1*v2*log(2)))/(2*q(k)*(v3-itr)^2)-2+lambda(k)*v4;
    end
    part_itr_total = 4*(itr-1)+ sum(part_itr(k));
    itr = max(round(itr-t_x*part_itr_total),1);

    %% 迭代更新拉格朗日乘子 lambda
    for k=1:m
        Grad_f(k) = p(k)*z*T*inv_pos(q(k)*itr*b*log(1+p(k)*hk/(I+bkn0))/log(2))+T*kapa*D*c*f(k)^2 - erequirement(k);
    end
    lambda = max( lambda + t_lambda * Grad_f' , zeros( m,1) );

    for k=1:m
    Term1(k) = 2^(z/(b*q(k)*(T_max(k)-itr*D*c/f(k)))-1);
    end
    obj = 2*(itr-1)^2;
    for k=1:m
        obj = obj + itr^2 * (exp((I+bkn0)*Term1(k)*inv_pos(p(k))/(distance(k)^(-gamma)))/q(k)-1);
    end
    
    % 拉格朗日函数
    L = obj+sum(Grad_f * lambda);
    if isinf(L)
        break;
        display('inf了')
    end
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
itr

% %% 查看收敛性
% figure(1)
% plot( 2:length(L_v),L_v(2:length(L_v)), '-ro', 'linewidth',2)
% grid on
% xlabel('Iteration round');
% ylabel('Lagrange function');
% 
% figure(2)
% plot( 2:length(L_v), obj_v(2:length(obj_v)) - L_v(2:length(L_v)) , '-b*', 'linewidth',2)
% grid on
% xlabel('Iteration round');
% ylabel('Objective - Lagrange function');
