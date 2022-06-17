function prob = prob_coll_n_beta(n,beta,i,T)
%prob_coll_n_beta - Description
%
% Syntax: prob = prob_coll_n_beta(n,beta,T,i)
%
% Long description
    prob = prob_n_beta(n,beta,i,T);
    if i ~= 1
        prob = prob .* prod(1 - prob_n_beta(n,beta,1:i-1,T));
    end
end