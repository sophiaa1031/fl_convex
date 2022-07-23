clc;clear;close

m=50; %Users' total number
iterationNumber=10;    %Maximun iteration times
hk=0.000004;    %Channel gain 
I=3*10^(-8);  %Interference
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

% f= 7e8+(1e9-7e8)*rand(1,m);
% distance = 1+12.*rand(1,m);
% v = 30.*rand(1,m);
cvx_begin
    variable p(1,m) nonnegative
    minimize(norm(exp((I+bkn0) .*inv_pos(p).*(2.^(z./(b.* q.*(T_max-itr*D*c./f))-1))./(distance.^(-gamma)))./q,1))
    subject to
       z*p(k)/(b* q(k)*log(1+p(k)*hk/(I+bkn0))/log(2))+itr*kapa*D*c*(f(k))^2<=erequirement;
       z/(b* q(k)*log(1+p(k)*hk/(I+bkn0))/log(2))+I*D*c/f(k)<=T_max;
cvx_end
