% clc;clear;close all
% 

q_min = ones(1,m)*0.001;% q���½�
q_max = ones(1,m)*1;% q���Ͻ�
get_gradient_q_left = zeros(1,m);
get_gradient_q_right = zeros(1,m);
delta_q = 0.01; % ���ַ�q������ֹͣ����

mu = ones(m,1)*0.05;   % �������ճ���mu��ʼ��
beta = 20;   % �������ճ���beta��ʼ��
real_ite=1; % ��ǰ��������
G_v = [];
mu_v = [];
beta_v = []; %�������̳�������ֵ
L_v = [];  % �������պ���
obj_v = [];
tol = 10^(-4); % �ж���������ֵ
t_x=0.01;      % ���߱��� x �ĸ��²���
t_mu=0.01;  % �������ճ��ӵĸ��²���
t_beta=0.01;      % �������ճ��ӵĸ��²���
max_iteration = 3e3;  % ����������

for i = 1 : max_iteration
    %% ���ַ��ҵ�����ֵ
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
    % ���������������ճ��� mu
    for k=1:m
        Grad_mu(k) = p(k)*z*inv_pos(q(k)*b*log(1+p(k)*hk/(I+bkn0))/log(2))+itr*kapa*D*c*f(k)^2 - erequirement;
    end
    mu = max( mu + t_mu * Grad_mu' , zeros( m,1) );
%     mu = mu + t_mu * Grad_mu';
    
    % ���������������ճ��� beta
    Grad_beta = sum(q)-N;
    beta = max( beta + t_beta * Grad_beta' , 0 );

%     for k=1:m
%     Term1(k) = 2^(z/(b*q(k)*(T_max-itr*D*c/f(k))))-1;
%     end
    obj = 0;
    for k=1:m
        obj = obj + exp((I+bkn0)*inv_pos(p(k))/(distance(k)^(-gamma))*(2^(z/(b*q(k)*(T_max-itr*D*c/f(k))))-1))/q(k);
    end
    

    % �������պ���
    L = obj+sum(Grad_mu * mu)+sum(Grad_beta * beta);
    if isinf(L)
        break;
        display('inf��')
    end
    % ��¼��ʷ���½��
    G_v = [G_v Grad_mu];%���ӵ��ݶ�
    mu_v = [mu_v mu];%���ӵ�ֵ
    beta_v = [beta_v beta];%���ӵ�ֵ
    obj_v = [obj_v obj]; % obj_v append
    L_v = [L_v L]; % L_v append
    %% ��ֹ�����ж�
    if (real_ite>2) && ( abs( obj - L ) < tol )  % k == 5000 %k>3&&(f_v(k)-f_v(k-1)<tol)
        break;
    end
    real_ite=real_ite+1;
end

%
% if i == max_iteration
%     disp(['����������������']);
% end
% disp(['��������Ϊ��',num2str(real_ite),'��']);
% disp(['����ֵ��Ϊ��']);
% q
% 
% % �鿴������
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
% %% p��f��d��ϵ
% figure(3)
% qfd=[(q*100)',(f/1e7)',distance'];
% bar(1:10,qfd(1:10,:))
% legend("q*10^(-2)","f*10^7","distance")
%      
