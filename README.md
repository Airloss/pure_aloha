# PureALOHA_CRP
A research about pure ALOHA with CRP
```MATLAB
clear

T = 100;
N = 1e5;

lambda = 0.3;

cnt_list = zeros(N,1);
s3_list = zeros(N,1);

for idx = 1:N
    pkt_list = exprnd(1/lambda,500,1);
    pkt_list = cumsum(pkt_list,1);
    conut = sum(pkt_list < 100);
    s3_list(idx) = pkt_list(3,1);
    cnt_list(conut,1) = cnt_list(conut,1) + 1;
end
cnt_list = cnt_list ./ N;
s3_list = sortrows(s3_list,1);
t = 0:(ceil(s3_list(end)) / N):ceil(s3_list(end));
s3_cdf = zeros(length(t),1);
for tdx = 1:length(t)
    cnt = sum(s3_list < t(tdx));
    s3_cdf(tdx) = cnt / length(t);
end

figure
plot(1:60,cnt_list(1:60,1),'--',1:60,poisspdf(1:60,lambda*T),'LineWidth',1)
legend('Simulation','Poisson','Interpreter','latex','Location','northwest','FontSize',14.4)
xlabel('Packets number','Interpreter','latex','FontSize',17.6)
ylabel('Probability','Interpreter','latex','FontSize',17.6)
title('Simulation of Poisson Process in $\lambda=0.3$','Interpreter','latex','FontSize',17.6)

figure
plot(t,s3_cdf,'--',t,F_s3(lambda .* t),'LineWidth',1)
legend('Simulation','CDF','Interpreter','latex','Location','northwest','FontSize',14.4)
xlabel('t','Interpreter','latex','FontSize',17.6)
ylabel('CDF Probability','Interpreter','latex','FontSize',17.6)
title('CDF of $S_3$','Interpreter','latex','FontSize',17.6)

function cdf_s3 = F_s3(lambda_t)
    cdf_s3 = 1 - exp(-lambda_t) .* (1 + lambda_t + lambda_t .^ 2 ./ 2);
end
```
