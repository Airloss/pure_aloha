function prob_j = prob_n_beta(n,beta,j,T)
%prob_n_beta - Description
%
% Syntax: prob_j = prob_n_beta(n,beta,j,T)
%
% Long description
    prob_j = exp(-(n-j) * beta * T);
end