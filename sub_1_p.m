clc;clear;close

m=50; %Users' total number
iterationNumber=10;    %Maximun iteration times
hk=0.000004;    %Channel gain 
I=3.7*10^(-8);  %Interference
bkn0=10^(-14);  %Noise = B * N0
b= 1e6; %bandwidth
gamma=3.76;
p0=0.01;  %uplink power
itr = 2;%local iterations
q = rand(1,m);%RB assignment
Dia = 15; %diameter of the BS
z = 28.1; %model size
T_max = 35*ones(1,m); %maximum latency in one round
erequirement=0.3*ones(1,m); %% energy requirement
kapa = 1e-28;%processor coefficient
D = 128;%sample size
c = 1.68e7;%cpu bits

f=normrnd(10^9,10^8,1,m);  %CPU frequency
distance = normrnd(7,2,1,m); %distance to the BS
v = normrnd(13,5,1,m); %velocity


p = ones(1,m)*0.1;
lambda = ones(m,1)*0.05;   %�������ճ��ӳ�ʼ�� 
real_ite=1;
G_v = [];  lambda_v = [];  L_v = [];  obj_v = [];
tol = 10^(-4);
t_lambda=0.01;  % �������ճ��ӵĸ��²���
t_x=0.1;      % ���߱��� x �ĸ��²���
max_iteration = 3e4;

theta = 1;
for i = 1 : max_iteration
    %% �������¹���P
    for k=1:m
        H_1 = (I+bkn0)*(2^(z/(b*q(k)*(T_max(k)-itr*D*c/f(k)))-1))/...
            (distance(k)^(-gamma));
        H_2 = z/q(k)/b;
        H_3 = hk/(I+bkn0);
        H_4 = -H_2*H_3*p(k)/log(2)/(1+H_3*p(k))/(log(1+H_3*p(k))/log(2))^2+...
            H_2/(log(1+H_3*p(k))/log(2));
        part_p = -H_1*theta*exp(H_1/p(k))+lambda(k)*H_4;
        
        p(k) = p(k)-t_x*part_p;
    end
    %% ���������������ճ��� lambda
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

obj = obj*theta;
L = obj+sum(Grad_f * lambda);
    %% ��¼��ʷ���½��
    G_v = [G_v Grad_f];%���ӵ��ݶ�
    lambda_v = [lambda_v lambda];%���ӵ�ֵ
    obj_v = [obj_v obj]; 
    L_v = [L_v L];
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

figure(1)
plot( 2:length(L_v),L_v(2:length(L_v)), '-ro', 'linewidth',2)

figure(2)
plot( 2:length(L_v), obj_v(2:length(obj_v)) - L_v(2:length(L_v)) , '-b*', 'linewidth',2)
