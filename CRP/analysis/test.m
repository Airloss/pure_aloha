clear

crp_l = zeros(100,1);
crp_s = crp_l;

crp_l(2) = crp_l_i(crp_l,0.25,2,1);
crp_s(2) = crp_s_i(crp_s,0.25,2,1);
crp_s(3) = crp_s_i(crp_s,0.25,3,1);
crp_l(3) = crp_l_i(crp_l,0.25,3,1);
crp_l(4) = crp_l_i(crp_l,0.25,4,1);
a = 1 - prob_n_beta(10,0.25,1:3,1);
b = prod(a);