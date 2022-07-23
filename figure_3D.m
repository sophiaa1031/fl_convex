%%
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

% f=normrnd(10^9,10^8,1,m);  %CPU frequency
% distance = normrnd(7,2,1,m); %distance to the BS
% v = normrnd(13,5,1,m); %velocity

% f= 1e9+2e9*rand(1,m);
% distance = 2+48*rand(1,m);
f = [2477280583.99080,2171974071.65295,1493469051.97195,2332832434.63894,1166965627.20525,2251919570.34317,2321889115.89469,2459503710.63444,2781504232.65064,2964606445.76721];
distance = [38.9133960961230,29.9094314180191,46.5590269910810,29.8443375564052,2.81518104018854,7.80125941273079,43.4101144975841,25.2462325381809,42.5530723796606,12.0514440330049];

epsilon1 = ones(1,m)*0.0001; %p的收敛条件
epsilon2 = ones(1,m)*0.0001; %q的收敛条件

g=4;
% T_max_loss = (0.6:0.01:0.9); %maximum latency in one round
T_max  = 0.6;
erequirement_loss=(0.0025:0.001:0.0055); %% energy requirement
erequirement = 0.004;
loss_function = zeros(4,g);

p_min = 0.001;
p_max = 0.01;
q_min = 0.004;
q_max = 1;
% p = rand(1,m)*(p_max-p_min)+p_min;
% q = rand(1,m)*(q_max-q_min)+q_min; 
p = [0.00562981110855134,0.00895852920814260,0.00629223449777648,0.00239277113790440,0.00279876540571707,0.00466259353425016,0.00773835146394122,0.00843025434207540,0.00810966726950078,0.00386671820859093];
q = [0.535927870861244,0.0935908760554987,0.115258921216431,0.139747378742545,0.679937695580988,0.497196311013302,0.192951564393510,0.497025801690260,0.151017789088782,0.0587542503185634]; 
itr = 5;%local iterations初始化
%%
for s = 1:4
    p = rand(1,m)*(p_max-p_min)+p_min;
    q = rand(1,m)*(q_max-q_min)+q_min;    
    for j = 1:g
%         if s ==2 %只优化p
%             p = ones(1,m)*0.1; % 功率p初始化
%             q = ones(1,m)*0.2;
%         elseif s ==3 %只优化q
%             p = ones(1,m)*0.1;
%             q = ones(1,m)*0.2;
%         elseif s ==4 %同时优化p,q
%             p = ones(1,m)*0.1;
%             q = ones(1,m)*0.2;%RB assignment初始化
%         end
        itr = 5;%local iterations初始化

        p_int = p; % 功率p初始化
        q_int = q;%RB assignment初始化
        itr_int = itr;%local iterations初始化
    
        p_last = zeros(1,m);
        q_last = zeros(1,m);
        round_total = 0;

        loss_temp = zeros(1,m);
    
    %     T_max = T_max_loss(j);
        T_max = 5;
        erequirement = erequirement_loss(j);
    
        %%
        if s > 1
%             if s >2
%                sub_2_q;%更新q
%             end
%         while(sum(abs(p_last-p)>=epsilon1)>0 && sum(abs(q_last-q)>=epsilon2)>0)
        while(round_total<5)
            
            q_last =q;
            if s ~= 2
            sub_2_q;%更新q
            end
            
            p_last = p;
            if s ~= 3
            sub_1_p; %更新p
%               p = p_int;
            end
        
%             sub_3_I;
    
            round_total = round_total+1;
        end
        end
        %%
        for k = 1:m
            loss_temp(k) =exp((I+bkn0)*inv_pos(p(k))/(distance(k)^(-gamma))*(2^(z/(b* q(k)*(T_max-itr*D*c/f(k))))-1));
        end
        %%
        loss_function(s,j) = sum(loss_temp ./ q);
        display(['s= ',num2str(s),'Round ',num2str(j),': loss function is ',num2str(loss_function(s,j))]);
        p
        q
    end 
end

 
figure;
plot(erequirement_loss,loss_function(1,:),'b-o',erequirement_loss,loss_function(2,:),'r-o',erequirement_loss,loss_function(3,:),'g-o',erequirement_loss,loss_function(4,:),'k-o','Linewidth',2,'MarkerSize',8);
grid on
% xlabel('maximum latency');
xlabel('maximum power');
ylabel('loss function upper bound');
legend('random','optimal p','optimal q','optimal p & q')