function coll_i = expect_coll_n_beta(n,beta,i,T)
%expect_n_beta - Description
%
% Syntax: coll_i = expect_n_beta(n,beta,i)
%
% Long description
    % for j = 1:i-1
    %     coll_i = coll_i + expect_intvl_n_beta(n,beta,j,T);
    % end
    coll_ = expect_intvl_n_beta(n,beta,1:i-1,T);
    coll_i = T + sum(coll_);
end

function intvl_j = expect_intvl_n_beta(n,beta,j,T)
%expect_intvl_n_beta - Description
%
% Syntax: intvl_j = expect_intvl_n_beta(n,beta,j,T)
%
% Long description
    intvl_j = (1 - (1 + (n-j)*beta*T) .* exp(-(n-j)*beta*T)) ./ ((1 - exp(-(n-j)*beta*T)) .* (n-j) * beta);
end