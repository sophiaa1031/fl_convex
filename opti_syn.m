clear_clc;

N = 5;T = N;H = 50;
Es=300;Eu = 3.6e9;B = 1e7;rou = -50;
sigma2 = -110;gamma = 1e-27;fc = 1.5e9;
vmax = 20;Ph = 10;M = 10;phi = rou/sigma2;
 
% dataSize
nodes_data = rand(1,N)*(1e5)*2;
% cycles
nodes_cyc = rand(1,N)*(1e10);
w = rand(N,2)*500;

q_temp=w;
% q_temp=rand(N,2)*500;
up_temp = ones(N,T);
a_temp=diag(ones(1,N)); %N½×µ¥Î»¾ØÕó
% a_temp=a_temp(randperm(N),:); %Ëæ¼´ÅÅÁÐ¾ØÕóµÄÐÐ

aoi_val=[];

aoi_p_last=1e6;
aoi_p=1e4;

aoi_last=1e6;
aoi_temp=1e4;

 while(aoi_last-aoi_temp>=1e-3)
    aoi_last = aoi_temp;
    opti_a_yalmip;
    aoi_val = [aoi_val aoi_a];
    
    aoi_p_last=1e6;
    aoi_p=1e4;
    while(aoi_p_last-aoi_p>1e-2)
        aoi_p_last = aoi_p;
        opti_p_dual_asend;
    end
    aoi_val = [aoi_val aoi_p];
    
    opti_q;
    aoi_val = [aoi_val aoi_q];
    
    aoi_temp = aoi_val(end);
end

