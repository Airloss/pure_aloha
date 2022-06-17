clear

%T = [(0.01:0.01:0.09),(0.1:0.05:1),(1.1:0.1:2)];
T = 0.01:0.01:1;

tic
for tdx = 1:length(T)
    num = 1000;
    beta_0 = 1 / (2 * 2 * T(tdx));
    
    crp_l = zeros(num,1);
    crp_s = crp_l;

    for idx = 2:num
        crp_l(idx) = crp_l_i(crp_l,beta_0,idx,T(tdx));
        crp_s(idx) = crp_s_i(crp_s,beta_0,idx,T(tdx));
    end
    fname1 = sprintf('data/finite/CRP_L_%.2f_%d.mat',T(tdx),num);
    fname2 = sprintf('data/finite/CRP_S_%.2f_%d.mat',T(tdx),num);
    save(fname1,'crp_l');
    save(fname2,'crp_s');
end
toc

function l_i = crp_l_i(l,beta,i,T)
% Compute crp(k)'s average duration
    switch i
    case 1
        disp ERROR_I
    case 2
        l_i = (1/(i*beta) + prob_coll_n_beta(i,beta,1,T) * (2*T+l(i-1)) + prob_coll_n_beta(i,beta,i,T) * expect_coll_n_beta(i,beta,i,T)) / (1 - prob_coll_n_beta(i,beta,i,T));
    otherwise
        prob_coll = zeros(i-2,1);
        expect_coll = prob_coll;
        for jdx = 2:i-1
            prob_coll(jdx-1) = prob_coll_n_beta(i,beta,jdx,T);
            expect_coll(jdx-1) = expect_coll_n_beta(i,beta,jdx,T);
        end
        sum_ = sum(prob_coll .* (expect_coll + l(2:i-1)));
        l_i = (1/(i*beta) + prob_coll_n_beta(i,beta,1,T) * (2*T+l(i-1)) + prob_coll_n_beta(i,beta,i,T) * expect_coll_n_beta(i,beta,i,T) + sum_) / (1 - prob_coll_n_beta(i,beta,i,T));
    end
end

function s_i = crp_s_i(s,beta,i,T)
% Compute crp(k)'s average success
    switch i
    case 1
        disp ERROR_I
    case 2
        s_i = 2;
    otherwise
        prob_coll = zeros(i-2,1);
        for jdx = 2:i-1
            prob_coll(jdx-1) = prob_coll_n_beta(i,beta,jdx,T);
        end
        sum_ = sum(prob_coll .* s(2:i-1));
        s_i = (prob_coll_n_beta(i,beta,1,T) * (1+s(i-1)) + sum_) / (1-prob_coll_n_beta(i,beta,i,T));
    end
end