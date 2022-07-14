clc;clear;close all

m=50; %Users' total number
N = m*0.5; %the number of resource blocks
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
T_max = 35; %maximum latency in one round
erequirement=0.7*ones(1,m); %% energy requirement
kapa = 1e-28;%processor coefficient
D = 128;%sample size
c = 1.68e7;%cpu bits

% f=normrnd(10^9,10^8,1,m);  %CPU frequency
% distance = normrnd(7,2,1,m); %distance to the BS
% v = normrnd(13,5,1,m); %velocity

f= 7e8+(1e9-7e8)*rand(1,m);
distance = 1+12.*rand(1,m);
v = 30.*rand(1,m);
p = rand(1,m); % ����p��ʼ��

q_min = ones(1,m)*0.001;% q���½�
q_max = ones(1,m)*1;% q���Ͻ�
q = ones(1,m)*1;% ���߱�����Դ����q��ʼ��
get_gradient_q_left = zeros(1,m);
get_gradient_q_right = zeros(1,m);
delta_q = 0.001; % ���ַ�q������ֹͣ����

mu = ones(m,1)*0.1;   % �������ճ���mu��ʼ��
beta = 10;   % �������ճ���beta��ʼ��
real_ite=1; % ��ǰ��������
G_v = [];
mu_v = [];
L_v = [];  % �������պ���
obj_v = [];
tol = 10^(-4); % �ж���������ֵ
t_mu=0.01;  % �������ճ��ӵĸ��²���
t_beta=0.01;      % �������ճ��ӵĸ��²���
max_iteration = 3e4;  % ����������

for i = 1 : max_iteration
    %% ���ַ��ҵ�����ֵ
    q_left = q_min;
    q_right = q_max;
    for k=1:m
      get_gradient_q_left(k) = get_gradient_q(q_left(k),f(k),distance(k),p(k),mu(k),beta);
      get_gradient_q_right(k) = get_gradient_q(q_right(k),f(k),distance(k),p(k),mu(k),beta);
      if get_gradient_q_left(k)<=0 && get_gradient_q_right(k)<=0
        q(k) = q_right(k);
      else
        if get_gradient_q_left(k)>=0 && get_gradient_q_right(k) >=0
          q(k) = q_left(k);
        else
          while q_right(k)-q_left(k)>=delta_q
            q(k) = (q_left(k)+q_right(k))/2;
            l_qk = get_gradient_q(q(k),f(k),distance(k),p(k),mu(k),beta);
            if (l_qk==0)
              break
            elseif (l_qk<0)
              q_left(k) = q(k);
            else
              q_right(k) = q(k);
            end
          end
        end
      end
    end
    % ���������������ճ��� mu
    for k=1:m
        Grad_mu(k) = p(k)*z*inv_pos(q(k)*b*log(1+p(k)*hk/(I+bkn0))/log(2))+itr*kapa*D*c*f(k)^2 - erequirement(k);
    end
    mu = max( mu + t_mu * Grad_mu' , zeros( m,1) );
    
    % ���������������ճ��� beta
    Grad_beta = sum(q)-N;
    beta = max( beta + t_beta * Grad_beta' , 0 );

    for k=1:m
    Term1(k) = 2^(z/(b*q(k)*(T_max-itr*D*c/f(k)))-1);
    end
    obj = 0;
    for k=1:m
        obj = obj + exp((I+bkn0)*Term1(k)*inv_pos(p(k))/(distance(k)^(-gamma)))/q(k);
    end
    
    % �������պ���
    L = obj+sum(Grad_mu * mu)+sum(Grad_beta * beta);
    % ��¼��ʷ���½��
    G_v = [G_v Grad_mu];%���ӵ��ݶ�
    mu_v = [mu_v mu];%���ӵ�ֵ
    obj_v = [obj_v obj]; % obj_v append
    L_v = [L_v L]; % L_v append
    %% ��ֹ�����ж�
    if (real_ite>2) && ( abs( obj - L ) < tol )  % k == 5000 %k>3&&(f_v(k)-f_v(k-1)<tol)
        break;
    end
    real_ite=real_ite+1;
end


if i == max_iteration
    disp(['����������������']);
end
disp(['��������Ϊ��',num2str(real_ite),'��']);
disp(['����ֵ��Ϊ��']);
q

%% �鿴������
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
