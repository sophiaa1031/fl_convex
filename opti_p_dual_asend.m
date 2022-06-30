
lambda = ones( N+1 , 1 )*5;   %拉格朗日乘子初始化 
k=1;
G_v = [];  lambda_v = [];  x_v = [];  L_v = [];  f_v = [];
tol = 10^(-4);
t_lambda=0.01;  % 拉格朗日乘子的更新步长
t_x=0.1;      % 决策变量 x 的更新步长
max_iteration = 3e4;

for i = 1 : max_iteration
    %% 迭代更新功率P
    for n=1:N
        for t=1:T
            Dist_nt_2 = (H^2+(q_temp(t,:)-w(n,:))*(q_temp(t,:)-w(n,:))');
            Term1 = -a_temp(n,t)*nodes_data(n)* up_temp(n,t)*lambda(n)*phi*log(2);
            Term2 = -a_temp(n,t)*nodes_data(n)*((1+N-t)/N+Ph*lambda(N+1))*phi*log(2);
            Term3 = B*Dist_nt_2*(1+ up_temp(n,t)*phi/Dist_nt_2)*(log(1+ up_temp(n,t)*phi/Dist_nt_2))^2;

            up_temp(n,t)= up_temp(n,t)-t_x*((Term1+Term2)/Term3+...
                        a_temp(n,t)*nodes_data(n)*lambda(n)*log(2)/(B*log(1+up_temp(n,t)*phi/Dist_nt_2)));
%             up_temp(n,t) = up_temp(n,t)-t_x*((-a_temp(n,t)*nodes_data(n)*lambda(n)*phi*log(2)+Term2)/Term3);
        end
    end
    
    %% 迭代更新拉格朗日乘子 lambda 
for n=1:N
    for t=1:T
        rate(n,t)=B*log(1+up_temp(n,t)*phi/( H*H + (q_temp(t,:)-w(n,:))*(q_temp(t,:)-w(n,:))' ))/log(2);
    end
end
 
e_f_p=0;
for t=1:(T-1)
    e_f_p = e_f_p+0.5*M*vmax*norm(q_temp(t+1,:)-q_temp(t,:));
end
e_f_p = e_f_p + 0.5*M*vmax*norm(q_temp(1,:)-[0 0]);
e_c_p=0;
for n=1:N
    e_c_p= e_c_p+gamma*fc*fc*nodes_cyc(n);
end

part_3 = 0;
for n=1:N
    for t=1:T
        part_3=part_3+Ph*a_temp(n,t)*(nodes_data(n)/rate(n,t)+nodes_cyc(n)/fc);
    end
end
part_3=part_3+e_f_p+e_c_p-Eu;

part_2_right = sum(diag(nodes_data)*a_temp.*up_temp./rate,2)-Es;

Grad_f = part_2_right;
Grad_f(N+1) = part_3;
lambda = max(   lambda + t_lambda * Grad_f  ,  zeros( N+1,1)    );  
    
part_1 = 0;
for n=1:N
    for t=1:T
        part_1=part_1+(N-t+1)/N*a_temp(n,t)*nodes_data(n)/rate(n,t);
    end
end

f = part_1;
L = part_1+sum(lambda.*[part_2_right;part_3]);

    %% 记录历史更新结果
    G_v = [G_v Grad_f];
%     x_v = [x_v x];
    lambda_v = [lambda_v lambda];
    f_v = [f_v f]; 
    L_v = [L_v L];
    
    %% 终止条件判断
    % if (k>2) &&( abs(L - L_v(k-1)) < tol )  % k == 5000 %k>3&&(f_v(k)-f_v(k-1)<tol)
    if (k>2) &&( abs( f - L ) < tol )  % k == 5000 %k>3&&(f_v(k)-f_v(k-1)<tol)
        break;
    end
    k=k+1;
end


if i == max_iteration
    disp(['超过最大迭代次数！']);
end
disp(['迭代次数为：',num2str(k),'次']);
disp(['最优值点为：']);
up_temp

% figure(1)
% plot( 2:length(L_v),L_v(2:length(L_v)), '-ro', 'linewidth',2)
% 
% figure(2)
% plot( 2:length(L_v), f_v(2:length(f_v)) - L_v(2:length(L_v)) , '-b*', 'linewidth',2)

aoi_part_1_p = 0;
for t=1:T-1
    aoi_part_1_p = aoi_part_1_p+(N-t)/N*norm(q_temp(t+1,:)-q_temp(t,:))/vmax;
end
aoi_part_1_p = aoi_part_1_p + norm(q_temp(1,:)-[0 0])/vmax;

aoi_part_2_2_p = 0;
for t=1:T
    for n=1:n
        aoi_part_2_2_p = aoi_part_2_2_p+(N-t+1)/N*a_temp(n,t)*nodes_cyc(n)/fc;
    end
end
aoi_p = aoi_part_1_p + aoi_part_2_2_p + f;
% aoi_p = aoi_part_1_p + f;
disp(['最优值为：',num2str(aoi_p),'']);