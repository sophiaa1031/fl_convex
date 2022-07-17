figure
q = [0.11:0.01:1];
y = zeros(1,90);
for i =1:90
    y(i) = get_gradient_q(q(i),5*10^8,15,0.5,0.1,10);
end
plot( q, y, '-b*', 'linewidth',2)
grid on
xlabel('Iteration round');
ylabel('Objective - Lagrange function');