%% clc;clear;close

lambda = ones(m,1)*0.5;   % 拉格朗日乘子lambda初始化 
mu = ones(m,1)*0.5;   % 拉格朗日乘子lambda初始化 
beta = ones(m,1)*0.5;   % 拉格朗日乘子lambda初始化 
real_ite=1; % 当前迭代次数
G_v = [];  
lambda_v = [];  
L_v = [];  % 拉格朗日函数
obj_v = []; 
tol = 10^(-5); % 判断收敛的阈值
t_lambda=0.01;  % 拉格朗日乘子的更新步长
t_mu=0.01;  % 拉格朗日乘子的更新步长
t_beta=0.01;  % 拉格朗日乘子的更新步长
t_x=0.01;      % 决策变量 x 的更新步长
max_iteration = 3e4;  % 最大迭代次数

for i = 1 : max_iteration
    %% 迭代更新功率P
    for k=1:m
        H1 = (I+bkn0)*(2^(z/(b*q(k)*(T_max-itr*D*c/f(k))))-1)/(distance(k)^(-gamma));
        H2 = z/(b*q(k));
        H3 = hk/(I+bkn0);
        H4 = -H2*H3*p(k)/(1+H3*p(k))/(log2(1+H3*p(k)))^2+H2/(log2(1+H3*p(k)));
        part_p = -H1*exp(H1/p(k))/(q(k)*(p(k))^2)+lambda(k)*H4;
        
        p(k) = max(p(k)-t_x*part_p,p_min)
    end
    %% 迭代更新拉格朗日乘子 lambda
    for k=1:m
        Grad_lambda(k) = p(k)*z/(q(k)*b*log2(1+p(k)*hk/(I+bkn0)))+itr*kapa*D*c*f(k)^2 - erequirement;
    end
    lambda = max( lambda + t_lambda * Grad_lambda' , zeros( m,1) );
    %% 迭代更新拉格朗日乘子 mu
    for k=1:m
        Grad_mu(k) = p_min-p(k);
    end
    mu = max( mu + t_mu * Grad_mu' , zeros( m,1) );
    %% 迭代更新拉格朗日乘子 beta
    for k=1:m
        Grad_beta(k) = p(k)-p_max;
    end
    beta = max( beta + t_beta * Grad_beta' , zeros( m,1) );
%     for k=1:m
%     Term1(k) = 2^(z/(b*q(k)*(T_max-itr*D*c/f(k)))-1);
%     end
    %%
    obj = 0;
    for k=1:m
        obj = obj + exp((I+bkn0)/p(k)/(distance(k)^(-gamma))*(2^(z/(b*q(k)*(T_max-itr*D*c/f(k))))-1))/q(k);
    end
    %%
    % 拉格朗日函数
    L = obj+sum(Grad_lambda * lambda)+sum(Grad_mu * mu)+sum(Grad_beta * beta);
    if isinf(L)
        break;
        display('inf了')
    end
    % 记录历史更新结果
%     G_v = [G_v Grad_f];%乘子的梯度
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
obj

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
pfd=[(p*100)',(f/1e8)',distance'];
bar(1:10,pfd(1:10,:))
legend("p*10^(-2)","f*10^8","distance")