% clc;clear;close all
% 

q_min = ones(1,m)*0.001;% q的下界
q_max = ones(1,m)*1;% q的上界
get_gradient_q_left = zeros(1,m);
get_gradient_q_right = zeros(1,m);
delta_q = 0.01; % 二分法q的搜索停止条件

mu = ones(m,1)*0.05;   % 拉格朗日乘子mu初始化
beta = 20;   % 拉格朗日乘子beta初始化
real_ite=1; % 当前迭代次数
G_v = [];
mu_v = [];
beta_v = []; %迭代过程乘子所有值
L_v = [];  % 拉格朗日函数
obj_v = [];
tol = 10^(-4); % 判断收敛的阈值
t_x=0.01;      % 决策变量 x 的更新步长
t_mu=0.01;  % 拉格朗日乘子的更新步长
t_beta=0.01;      % 拉格朗日乘子的更新步长
max_iteration = 3e3;  % 最大迭代次数

for i = 1 : max_iteration
    %% 二分法找到最优值
%     q_left = q_min;
%     q_right = q_max;
    for k=1:m
%       get_gradient_q_left(k) = get_gradient_q(q_left(k),f(k),distance(k),p(k),mu(k),beta);
%       get_gradient_q_right(k) = get_gradient_q(q_right(k),f(k),distance(k),p(k),mu(k),beta);
%       if get_gradient_q_left(k)<=0 && get_gradient_q_right(k)<=0
%         q(k) = q_right(k);
%       else
%         if get_gradient_q_left(k)>=0 && get_gradient_q_right(k) >=0
%           q(k) = q_left(k);
%         else
%           while q_right(k)-q_left(k)>=delta_q
%             q(k) = (q_left(k)+q_right(k))/2;
%             l_qk = get_gradient_q(q(k),f(k),distance(k),p(k),mu(k),beta);
%             if (l_qk==0)
%               break
%             elseif (l_qk<0)
%               q_left(k) = q(k);
%             else
%               q_right(k) = q(k);
%             end
%           end
%         end
%       end
        V1 = (I+bkn0)/(p(k)*distance(k)^(-gamma));
        V2 = 2^(z/(b*(T_max-itr*D*c/f(k))));
        V3 = p(k)*z/(q(k)^2*b*log2(1+p(k)*hk/(I+bkn0)));
        l_q = -exp((2*(V2/q(k)-1))*V1)/q(k)^2-V1*V2*log(2)*2^(V2/q(k)-1)*exp((2*V2/q(k)-1)*V1)/q(k)^3-mu(k)*V3+beta;
        q(k) = max(min(q(k)-t_x*l_q,1),0);
    end
    % 迭代更新拉格朗日乘子 mu
    for k=1:m
        Grad_mu(k) = p(k)*z*inv_pos(q(k)*b*log(1+p(k)*hk/(I+bkn0))/log(2))+itr*kapa*D*c*f(k)^2 - erequirement;
    end
    mu = max( mu + t_mu * Grad_mu' , zeros( m,1) );
%     mu = mu + t_mu * Grad_mu';
    
    % 迭代更新拉格朗日乘子 beta
    Grad_beta = sum(q)-N;
    beta = max( beta + t_beta * Grad_beta' , 0 );

%     for k=1:m
%     Term1(k) = 2^(z/(b*q(k)*(T_max-itr*D*c/f(k))))-1;
%     end
    obj = 0;
    for k=1:m
        obj = obj + exp((I+bkn0)*inv_pos(p(k))/(distance(k)^(-gamma))*(2^(z/(b*q(k)*(T_max-itr*D*c/f(k))))-1))/q(k);
    end
    

    % 拉格朗日函数
    L = obj+sum(Grad_mu * mu)+sum(Grad_beta * beta);
    if isinf(L)
        break;
        display('inf了')
    end
    % 记录历史更新结果
    G_v = [G_v Grad_mu];%乘子的梯度
    mu_v = [mu_v mu];%乘子的值
    beta_v = [beta_v beta];%乘子的值
    obj_v = [obj_v obj]; % obj_v append
    L_v = [L_v L]; % L_v append
    %% 终止条件判断
    if (real_ite>2) && ( abs( obj - L ) < tol )  % k == 5000 %k>3&&(f_v(k)-f_v(k-1)<tol)
        break;
    end
    real_ite=real_ite+1;
end

%
% if i == max_iteration
%     disp(['超过最大迭代次数！']);
% end
% disp(['迭代次数为：',num2str(real_ite),'次']);
% disp(['最优值点为：']);
% q
% 
% % 查看收敛性
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
% 
% %% p与f，d关系
% figure(3)
% qfd=[(q*100)',(f/1e7)',distance'];
% bar(1:10,qfd(1:10,:))
% legend("q*10^(-2)","f*10^7","distance")
%      
