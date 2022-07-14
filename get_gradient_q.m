
function [ l_q ] = get_gradient_p(q,f,distance,p,mu,beta)
    
m=30; %Users' total number
N = m*0.5; %the number of resource blocks
iterationNumber=10;    %Maximun iteration times
hk=0.000004;    %Channel gain
I=3.7*10^(-8);  %Interference
bkn0=10^(-14);  %Noise = B * N0
b= 1e6; %bandwidth
gamma=3.76;
itr = 2;%local iterations
z = 28.1; %model size
T_max = 35; %maximum latency in one round
erequirement=0.7*ones(1,m); %% energy requirement
D = 128;%sample size
c = 1.68e7;%cpu bits


    V_1 = (I+bkn0)/(p*distance^(-gamma));
    V_2 = 2^(z/(b*(T_max-itr*D*c/f)));
    V_3 = p*z/(q^2*b*log2(1+p*hk/(I+bkn0)));
    l_q = -exp((2*(V_2/q-1))*V_1)/q^2-V_1*V_2*log(2)*2^(V_2/q-1)*exp((2*V_2/q-1)*V_1)/q^3-mu*V_3+beta;
end

