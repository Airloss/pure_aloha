clear

THEATA = 0.99;
ENDTIME = 2e5;

lambda = 0.01:0.01:0.36;

thrpt_list = zeros(length(lambda),1);
dly_list = zeros(length(lambda),1);
scs_crp_list = zeros(length(lambda),1); % Record the mean of success packets in CRP

tic
parfor (ldx = 1:length(lambda),6)
    num = ceil(3 * lambda(ldx) * ENDTIME);
    ptr = 1;
    blg = 0;
    scs = 0;
    dly = 0;
    cnt = 0;
    cnt_crp = 0;
    min_t = 0;
    es_blg = 10;
    ac_blg = 10;
    scs_crp = 0;
    lambda_recur = 0;
    prev_end_t = 0;
    mu = 0.5 / es_blg;

    pkt_list = zeros(num,3);
    pkt_list(:,1) = cumsum(exprnd(1/lambda(ldx), num, 1));
    pkt_list(:,2) = pkt_list(:,1);

    while min_t < ENDTIME
        cnt = cnt + 1;

        if blg == 0
            pkt_list(ptr,1) = pkt_list(ptr,1) + exprnd(mu,1);
            pkt_list(ptr,3) = 1;    % stack in backlog list
            blg = 1;
        end

        min_t_temp = pkt_list(ptr,1) + 1;    % packet length equals 1
        new_pkt = sum(pkt_list(ptr+blg:end,1) < min_t_temp);
        if new_pkt > 0
            bof = exprnd(1/mu, new_pkt, 1);
            min_t_temp = min(min(pkt_list(ptr+blg:scs+blg+new_pkt,1)+bof)+1, min_t_temp);
            new_blg = sum(pkt_list(ptr+blg:scs+blg+new_pkt,1) < min_t_temp);
            pkt_list(ptr+blg:scs+blg+new_blg,1) = pkt_list(ptr+blg:scs+blg+new_blg,1) + bof(1:new_blg);
            pkt_list(ptr+blg:scs+blg+new_blg,3) = 1;
            blg = blg + new_blg;
            pkt_list(ptr:scs+blg,:) = sortrows(pkt_list(ptr:scs+blg,:),1);
        end
        idle_t = pkt_list(ptr,1) - min_t;
        if idle_t < 0
            disp FALSE_IDLE_MINUS
        end
        min_t = pkt_list(ptr,1) + 1;

        % cycle period
        sect = sum(pkt_list(ptr:scs+blg,1) < min_t);
        if sect == 1
            lambda_recur = THEATA * lambda_recur + (1-THEATA) / (min_t - prev_end_t);
            es_blg = max(es_blg * exp(-mu * idle_t) + lambda_recur,1);
            scs = scs + 1;
            dly = dly + pkt_list(ptr,1) - pkt_list(ptr,2) + 1;
            pkt_list(ptr,3) = -1;   % column 3 == -1 => scs
            ptr = ptr + 1;
            blg = blg - 1;
        else
            coll_start_t = pkt_list(ptr,1);
            lambda_recur = lambda_recur * THEATA;
            blg_end = sect - 1;

            % collision period
            while blg_end > 0 && min_t < ENDTIME
                min_t = pkt_list(ptr+blg_end,1) + 1;
                new_blg = sum(pkt_list(ptr+blg:end,1) < min_t);
                if new_blg > 0
                    pkt_list(ptr+blg:scs+blg+new_blg,1) = pkt_list(ptr+blg:scs+blg+new_blg,1) + exprnd(1/mu, new_blg, 1);
                    pkt_list(ptr+blg:scs+blg+new_blg,3) = 1;
                    blg = blg + new_blg;
                    pkt_list(ptr:scs+blg,:) = sortrows(pkt_list(ptr:scs+blg,:),1);
                end
                sect = sum(pkt_list(ptr:scs+blg,1) < min_t);
                if sect == blg_end + 1
                    break;
                else
                    blg_end = sect - 1;
                end
            end

            % CRP period
            cnt_crp = cnt_crp + 1;
            pkt_list(ptr:ptr+blg_end,1) = min_t;
            crp_start_t = min_t;
            crp_list = pkt_list(ptr:ptr+blg_end,:);
            if sum(crp_list(:,3) == 1) <= 1
                disp FALSE_COLL
            end
            crp_min_t = 0;
            crp_mu = 0.5 / 2;
            crp_ptr = 1;
            crp_blg = sum(crp_list(:,3)==1);
            crp_cnt = 0;
            crp_list(:,1) = crp_list(:,1) + exprnd(1/crp_mu,crp_blg,1);
            crp_list = sortrows(crp_list,1);
            while sum(crp_list(:,3) == 1)  > 0 && crp_min_t < ENDTIME
                if sum(crp_list(:,3) == 1) == 1
                    idx = find(crp_list(:,3) == 1);
                    if idx ~= crp_ptr
                        disp FASLSE_CRP_IDX
                    end
                    crp_list(idx,3) = -1;   % column 3 == -1 => scs
                    scs = scs + 1;
                    scs_crp = scs_crp + 1;
                    dly = dly + crp_list(idx,1) - crp_list(idx,2) + 1;
                    blg = blg - 1;
                    crp_ptr = crp_ptr + 1;
                    crp_min_t = crp_list(idx,1) + 1;
                else
                    crp_min_t = crp_list(crp_ptr,1) + 1;
                    crp_sect = crp_list(crp_ptr:end,1) < crp_min_t;
                    if sum(crp_sect) == 1
                        scs = scs + 1;
                        scs_crp = scs_crp + 1;
                        dly = dly + crp_list(crp_ptr,1) - crp_list(crp_ptr,2) + 1;
                        crp_list(crp_ptr,3) = -1;
                        crp_ptr = crp_ptr + 1;
                        crp_cnt = crp_cnt + 1;
                        blg = blg - 1;
                        if crp_cnt >= min(crp_blg - 1, 3)
                            break;
                        end
                    else
                        crp_sect_num = sum(crp_sect) - 1;
                        while crp_sect_num > 0 && crp_min_t < ENDTIME
                            crp_list(crp_ptr:crp_ptr+crp_sect_num-1,1) = crp_list(crp_ptr:crp_ptr+crp_sect_num-1,1) + 1 + exprnd(1/crp_mu,crp_sect_num,1);
                            crp_min_t = crp_list(crp_ptr+crp_sect_num,1) + 1;
                            crp_list(crp_ptr:end,:) = sortrows(crp_list(crp_ptr:end,:),1);
                            crp_sect_num = sum(crp_list(crp_ptr:end,1) < crp_min_t) - 1;
                        end
                        crp_list(crp_ptr,1) = crp_list(crp_ptr,1) + 1 + exprnd(1/crp_mu);
                        crp_list(crp_ptr:end,:) = sortrows(crp_list(crp_ptr:end,:),1);
%                        crp_cnt = 0;
                    end
                end
            end

            min_t = crp_min_t;
            coll_t = min_t - coll_start_t;
            es_blg = max(1 + es_blg * exp(-mu * idle_t) + lambda_recur * coll_t,1);
            pkt_list(ptr:ptr+blg_end,:) = crp_list;
            ptr = scs + 1;
            if sum(pkt_list(1:scs,3) > 0) > 0
                disp FALSE_CRP
            end
            old_blg = pkt_list(ptr:scs+blg,1) < min_t;
            if sum(old_blg > 0)
                pkt_list(ptr:scs+blg,1) = pkt_list(ptr:scs+blg,1) .* ~old_blg + (min_t + exprnd(1/mu,blg,1)) .* old_blg;
                pkt_list(ptr:scs+blg,:) = sortrows(pkt_list(ptr:scs+blg,:),1);
            end
            new_blg = sum(pkt_list(ptr+blg:end,1) < min_t);
            if new_blg > 0
                pkt_list(ptr+blg:scs+blg+new_blg,1) = min_t + exprnd(1/mu,new_blg,1);
                pkt_list(ptr+blg:scs+blg+new_blg,3) = 1;
                blg = blg + new_blg;
            end
            pkt_list(ptr:scs+blg,:) = sortrows(pkt_list(ptr:scs+blg,:),1);
            if sum(pkt_list(:,2) < min_t) - scs ~= blg
                disp FALSE_BLG1
            end
        end
        ac_blg = sum(pkt_list(:,2) < min_t) - scs;
        if ac_blg ~= blg
            disp FALSE_BLG
        end
        prev_end_t = min_t;
        mu = 0.5 / es_blg;
    end
    thrpt_list(ldx) = scs / min_t;
    dly_list(ldx) = dly / scs;
    scs_crp_list(ldx) = scs_crp / cnt_crp;
end
toc

figure
plot(lambda,thrpt_list,'LineWidth',1)
legend('Estimated','Location','northwest','Interpreter','latex','FontSize',14.4)
grid on
xlim([0 0.26])
xlabel('$\lambda$','Interpreter','latex','FontSize',17.6)
ylabel('Throughput (packet/sec)','Interpreter','latex','FontSize',17.6)
title('Pure ALOHA CRP','Interpreter','latex','FontSize',17.6)

figure
plot(lambda,dly_list,'LineWidth',1)
legend('Estimated','Location','northwest','Interpreter','latex','FontSize',14.4)
grid on
xlim([0 0.26])
ylim([0 300])
xlabel('$\lambda$','Interpreter','latex','FontSize',17.6)
ylabel('Delay (sec)','Interpreter','latex','FontSize',17.6)
title('Pure ALOHA CRP','Interpreter','latex','FontSize',17.6)
