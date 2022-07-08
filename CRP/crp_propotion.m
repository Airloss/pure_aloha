clear

USER = 3;
ENDTIME = 2e5;

lambda = 0.02:0.02:0.6;

crp_proportion = zeros(length(lambda),1);

for ldx = 1:length(lambda)
    ptr = ones(USER,1);
    scs = 0;
    dly = 0;
    trans_flag = 0;
    min_t = 0;
    crp_t = 0;

    num = ceil(1.5 * lambda(ldx) * ENDTIME);
    pkt_list = zeros(num,3,USER);
    trans_list = zeros(USER,3);
    for i1 = 1:USER
        pkt_list(:,1,i1) = cumsum(exprnd(1/lambda(ldx), num, 1));
        pkt_list(:,2,i1) = pkt_list(:,1,i1);
        pkt_list(:,3,i1) = i1;
        trans_list(i1,:) = pkt_list(1,:,i1);
    end
    trans_list = sortrows(trans_list,1);
       
    while min_t <= ENDTIME
        trans_flag = trans_list(1,3);
        if min_t < trans_list(1,1)
            min_t = trans_list(1,1) + 1;
        else
            min_t = min_t + 1;
            trans_list(:,1) = min_t;
        end
        if (trans_list(1,1) - trans_list(2,1)) >= 1
            crp_t = crp_t + 1;
            scs = scs + 1;
            dly = dly + min_t - trans_list(1,2);
            ptr(trans_flag) = ptr(trans_flag) + 1;
            trans_list(1,:) = pkt_list(ptr(trans_flag),:,trans_flag);
            trans_list = sortrows(trans_list,1);
        else
            num = sum(trans_list(:,1) <= min_t) - 1;
            if num < USER - 1
                while num > 0 && min_t <= ENDTIME
                    min_t = 
                end
            end
            
            crp_start_t = min_t;
            crp_list = trans_list;
            crp_mu = 0.36;
            crp_min_t = min_t;

            crp_list(:,1) = crp_list(:,1) + exprnd(1/crp_mu,2,1);
            crp_list = sortrows(crp_list,1);
            while sum(crp_list(:,3) >= 1) > 0 && crp_min_t < ENDTIME
                if sum(crp_list(:,3) >= 1) == 1
                    ic = find(crp_list(:,3) >= 1);
                    crp_list(ic,3) = - crp_list(ic,3);   % column 3 == -1 => scs
                    scs = scs + 1;
                    dly = dly + crp_list(ic,1) - crp_list(ic,2) + 1;
                    crp_min_t = crp_list(ic,1) + 1;
                else
                    crp_trans_idx = crp_list(:,3) >= 1;
                    crp_trans_list = crp_list(crp_trans_idx,:);
                    crp_min_t = crp_trans_list(1,1) + 1;
                    crp_sect = sum(crp_trans_list(:,1) < crp_min_t);
                    if crp_sect == 1
                        scs = scs + 1;
                        dly= dly + crp_trans_list(1,1) - crp_trans_list(1,2) + 1;
                        crp_trans_list(1,3) = - crp_trans_list(1,3);
                        crp_trans_list(2:end,1) = crp_trans_list(1,1) + 1;
                    else
                        crp_sect_num = crp_sect - 1;
                        while crp_sect_num > 0 && crp_min_t < ENDTIME
                            crp_min_t = crp_trans_list(crp_sect,1) + 1;
                            crp_sect = sum(crp_trans_list(:,1) < crp_min_t);
                            if crp_sect == crp_sect_num + 1
                                break;
                            else
                                crp_sect_num = crp_sect - 1;
                            end
                        end
                        crp_trans_list(1:crp_sect,1) = crp_min_t + exprnd(1/crp_mu,crp_sect,1);
                        crp_trans_list = sortrows(crp_trans_list,1);
                    end
                    crp_list(crp_trans_idx,:) = crp_trans_list;
                end
            end
            crp_t = crp_t + crp_min_t - crp_start_t;
            min_t = crp_min_t;
            crp_trans_list(:,3) = - crp_trans_list(:,3);
            ptr = ptr + 1;
            trans_list(1,:) = pkt_list(ptr(trans_flag),:,trans_flag);
            trans_list(2,:) = pkt_list(ptr(3 - trans_flag),:,3 - trans_flag);
            trans_list = sortrows(trans_list,1);
        end
    end
    crp_proportion(ldx) = crp_t / min_t;
end

figure
plot(lambda,crp_proportion,'LineWidth',1)
legend('CRP Propotion','Location','southeast','Interpreter','latex','FontSize',14.4)
grid on
xlabel('$\lambda$','Interpreter','latex','FontSize',17.6)
ylabel('Channel Efficiency','Interpreter','latex','FontSize',17.6)
% title('Pure ALOHA CRP','Interpreter','latex','FontSize',17.6)