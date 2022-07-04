clear

fname1 = sprintf('data/finite/CRP_L_%.2f_%d.mat',1,1000);
fname2 = sprintf('data/finite/CRP_S_%.2f_%d.mat',1,1000);
load(fname1,'crp_l');
load(fname2,'crp_s');
crp_l(1) = [];
crp_s(1) = [];

idx = 1:30;

figure
plot(idx,crp_l(idx),'LineWidth',1)
legend('Main Channel Throughput','Location','southeast','Interpreter','latex','FontSize',14.4)
grid on
% xlim([0 0.36])
xlabel('$k$','Interpreter','latex','FontSize',17.6)
ylabel('$L_{k}$','Interpreter','latex','FontSize',17.6)
title('CRP(k) Duration','Interpreter','latex','FontSize',17.6)