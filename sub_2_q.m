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

% f=normrnd(10^9,10^8,1,m);  %CPU frequency
% distance = normrnd(7,2,1,m); %distance to the BS
% v = normrnd(13,5,1,m); %velocity

f= 7e8+(1e9-7e8)*rand(1,m);
distance = 1+12.*rand(1,m);
v = 30.*rand(1,m);

q = ones(1,m)*0.1;% ���߱�����Դ����q��ʼ��
p = rand(1,m); % ����p��ʼ��

mu = ones(m,1)*0.05;   % �������ճ���mu��ʼ�� 
beta = ones(m,1)*0.05;   % �������ճ���beta��ʼ�� 
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
    
    %% ����д����
    
    
    %% ���������������ճ��� mu
    for k=1:m
        % �ݶ���Ҫ��
        Grad_f(k) = p(k)*z*inv_pos(q(k)*b*log(1+p(k)*hk/(I+bkn0))/log(2))+itr*kapa*D*c*f(k)^2 - erequirement(k);
    end
    mu = max( mu + t_mu * Grad_f' , zeros( m,1) );

    for k=1:m
    Term1(k) = 2^(z/(b*q(k)*(T_max(k)-itr*D*c/f(k)))-1);
    end
    obj = 0;
    for k=1:m
        obj = obj + exp((I+bkn0)*Term1(k)*inv_pos(p(k))/(distance(k)^(-gamma)))/q(k);
    end
    
    % �������պ���
    L = obj+sum(Grad_f * mu);
    % ��¼��ʷ���½��
    G_v = [G_v Grad_f];%���ӵ��ݶ�
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
p

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

%% p��f��d��ϵ
figure(3)
pfd=[(p*100)',(f/1e8)',distance']
bar(1:10,pfd(1:10,:))
legend("p*10^(-2)","f*10^8","distance")