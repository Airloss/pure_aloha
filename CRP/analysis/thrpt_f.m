%% Numerical search beta that satisfied the max throughput
% N*beta*T
clear

T = 1;
num = 200;
beta = 0.0001:0.0001:0.02;

thrpt = zeros(length(beta),1);
prob_cc = thrpt;

% load crp data
fname1 = sprintf('data/finite/CRP_L_%.2f_%d.mat',T,num);
fname2 = sprintf('data/finite/CRP_S_%.2f_%d.mat',T,num);
load(fname1,'crp_l');
load(fname2,'crp_s');
crp_l(1) = [];
crp_s(1) = [];

tic
parfor (bdx = 1:length(beta),6)
    prob_coll = zeros(num,1);
    expect_coll = prob_coll;
    prob_coll(1) = prob_n_beta(num,beta(bdx),1,T);  % success
    expect_coll(1) = T; % success
    % collision
    for idx = 2:num
        prob_coll(idx) = prob_coll_n_beta(num,beta(bdx),idx,T);
        expect_coll(idx) = expect_coll_n_beta(num,beta(bdx),idx,T);
    end
    expect_busy = prob_coll(1) * expect_coll(1) + sum(prob_coll(2:num) .* (expect_coll(2:num) + crp_l));
    thrpt(bdx) = (prob_coll(1) + sum(prob_coll(2:num) .* crp_s)) / (1 / (num * beta(bdx)) + expect_busy);
    prob_cc(bdx) = (1-prob_coll(1))^2;
end
toc

pt = find(thrpt == max(thrpt));
disp(beta(pt));
disp(thrpt(pt));
disp(beta(pt) * num * T);
disp(thrpt(pt) * T);

ftitle = sprintf('$N=%d$',num);
figure
plot(beta,thrpt,'LineWidth',1)
grid on
legend(ftitle,'Location','northeast','Interpreter','latex','FontSize',14.4)
xlabel('$\beta$ (packets/sec)','Interpreter','latex','FontSize',17.6)
ylabel('Throughput (packets/sec)','Interpreter','latex','FontSize',17.6)
title('$T=1$','Interpreter','latex','FontSize',17.6)