%% Numerical search GT that satisfied the max throughput
% GT = 0.64687, ST = 0.26032
clear

num = 1000;
%T = 0.01:0.01:1;
T = [(0.01:0.01:0.09),(0.1:0.05:1),(1.1:0.1:2)];
G = [(0.0001:0.0001:1),(1:0.1:3)];
thrpt = zeros(length(T),1);
gt = thrpt;

tic
for tdx = 1:length(T)
    fname1 = sprintf('data/finite/CRP_L_%.2f_%d.mat',T(tdx),num);
    fname2 = sprintf('data/finite/CRP_S_%.2f_%d.mat',T(tdx),num);
    load(fname1,'crp_l');
    load(fname2,'crp_s');
    crp_l(1) = [];
    crp_s(1) = [];
    thrpt_ = zeros(length(G),1);
    for gdx = 1:length(G)
        prob_coll = zeros(num,1);
        expect_coll = prob_coll;
        prob_coll(1) = exp(-G(gdx) * T(tdx));
        expect_coll(1) = T(tdx);
        prob_coll(2:num,1) = prob_coll_g(G(gdx), 2:num, T(tdx));
        expect_coll(2:num,1) = expect_coll_g(G(gdx), 2:num, T(tdx));
        expect_busy = prob_coll(1) * expect_coll(1) + sum(prob_coll(2:num) .* (expect_coll(2:num) + crp_l));
        thrpt_(gdx) = (prob_coll(1) + sum(prob_coll(2:num) .* crp_s)) / (1 / G(gdx) + expect_busy);
    end
    pt = find(thrpt_ == max(thrpt_));
    thrpt(tdx) = thrpt_(pt(1)) * T(tdx);
    gt(tdx) = G(pt(1)) * T(tdx);
end
toc

ftitle = sprintf('$N=%d$',num);
figure
plot(T,thrpt,T,gt,'LineWidth',1)
grid on
legend('$S_\infty T$','$GT$','Location','northeast','Interpreter','latex','FontSize',14.4)
xlabel('$T$ (sec)','Interpreter','latex','FontSize',17.6)
ylabel('$S_\infty T/GT$','Interpreter','latex','FontSize',17.6)
title(ftitle,'Interpreter','latex','FontSize',17.6)

function prob = prob_coll_g(G,i,T)
% Long description
    prob = exp(-G*T) .* (1 - exp(-G*T)).^(i-1);
end

function expect = expect_coll_g(G,i,T)
% Long description
    expect = T + (i-1) .* (1 - (1 + G*T) * exp(-G*T)) / ((1 - exp(-G*T)) * G);
end