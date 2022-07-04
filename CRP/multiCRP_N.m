clear

THEATA = 0.99;
ENDTIME = 1e5;
m = 2;

lambda = 0.3:0.1:0.8;

thrpt_list = zeros(length(lambda),1);
dly_list = zeros(length(lambda),1);
avg_coll = zeros(length(lambda),1);
avg_idle = zeros(length(lambda),1);

crp_scs_list = zeros(length(lambda),1); % Record the mean of success packets in CRP
crp_len = zeros(length(lambda),1);
crp_avg_coll = zeros(length(lambda),1);

crp_mu = 0.5 / 2;
lambda_recur = 0.2;

tic
parfor (ldx = 1:length(lambda),6)
    num = ceil(ENDTIME / (1/lambda(ldx)));
    ptr = 1;
    blg = 0;
    scs = 0;
    dly = 0;
    cnt = 0;
    min_t = 0;
    es_blg = 10;
    ac_blg = 10;
    coll_t = 0;
    idle_t = 0;
    cnt_coll = 0;
    crp_t = 0;
    crp_cnt = 0;
    crp_scs = 0;
    crp_coll = 0;
    prev_end_t = 0;
    mu = 0.6468 / es_blg;

    channel_flag = zeros(m,1);
    channel_mint = channel_flag;
    pkt_list = zeros(num,3);
    idx = 1;
    while pkt_list(idx,1) < ENDTIME && idx <= num
        sft = randi([1,ceil(lambda(ldx)*10)]);
        if idx + sft > num
            break;
        end
        pkt_list(idx:idx+sft-1,1) = pkt_list(idx,1) + exprnd(1/lambda(ldx),sft,1);
        pkt_list(idx + sft,1) = max(pkt_list(idx:idx+sft,1));
        idx = idx + sft;
    end
    pkt_list(idx:end,:) = [];
    pkt_list = sortrows(pkt_list,1);
    pkt_list(:,2) = pkt_list(:,1);

    while min_t < ENDTIME && ptr < num - 10
        cnt = cnt + 1;

        if blg == 0
            pkt_list(ptr,1) = pkt_list(ptr,1) + 0;
            pkt_list(ptr,3) = 1;    % stack in backlog list
            blg = 1;
        end

        min_t_temp = pkt_list(ptr,1) + 1;    % packet length equals 1
        new_pkt = sum(pkt_list(ptr+blg:end,1) < min_t_temp);
        if new_pkt > 0
            bof = 0;
            min_t_temp = min(min( ...
                pkt_list(ptr+blg:scs+blg+new_pkt,1)+bof)+1, min_t_temp);
            new_blg = sum(pkt_list(ptr+blg:scs+blg+new_pkt,1) < min_t_temp);
            pkt_list(ptr+blg:scs+blg+new_blg,1) = pkt_list( ...
                ptr+blg:scs+blg+new_blg,1) + bof;
            pkt_list(ptr+blg:scs+blg+new_blg,3) = 1;
            blg = blg + new_blg;
            pkt_list(ptr:scs+blg,:) = sortrows(pkt_list(ptr:scs+blg,:),1);
        end
        idle_t = idle_t + pkt_list(ptr,1) - min_t;
        if idle_t < 0
            disp FALSE_IDLE_MINUS
            % return
        end
        min_t = pkt_list(ptr,1) + 1;

        % cycle period
        sect = sum(pkt_list(ptr:scs+blg,1) < min_t);
        if sect == 1
            scs = scs + 1;
            dly = dly + pkt_list(ptr,1) - pkt_list(ptr,2) + 1;
            pkt_list(ptr,3) = -1;   % column 3 == -1 => scs
            ptr = ptr + 1;
            blg = blg - 1;
        else
            coll_start_t = pkt_list(ptr,1);
            blg_end = sect - 1;
            cnt_coll = cnt_coll + 1;
            % collision period
            while blg_end > 0 && min_t < ENDTIME
                min_t = pkt_list(ptr+blg_end,1) + 1;
                new_blg = sum(pkt_list(ptr+blg:end,1) < min_t);
                if new_blg > 0
                    pkt_list(ptr+blg:scs+blg+new_blg,1) = pkt_list( ...
                        ptr+blg:scs+blg+new_blg,1) + 0;
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
            coll_t = coll_t + min_t - coll_start_t;

            % CRP period in addtional channel
            flag = flag_status(channel_flag);
            if flag > 0
                channel_flag(flag) = 1;
                crp_cnt = crp_cnt + 1;
                pkt_list(ptr:ptr+blg_end,1) = min_t;
                crp_start_t = min_t;
                crp_list = pkt_list(ptr:ptr+blg_end,:);
                if sum(crp_list(:,3) == 1) <= 1
                    disp FALSE_COLL
                end
                crp_min_t = 0;
                crp_blg = sum(crp_list(:,3)==1);
                crp_list(:,1) = crp_list(:,1) + exprnd(1/crp_mu,crp_blg,1);
                crp_list = sortrows(crp_list,1);
                crp_coll = crp_coll + crp_blg;
                while sum(crp_list(:,3) == 1)  > 0 && crp_min_t < ENDTIME
                    if sum(crp_list(:,3) == 1) == 1
                        idx = find(crp_list(:,3) == 1);
                        crp_list(idx,3) = -1;   % column 3 == -1 => scs
                        scs = scs + 1;
                        crp_scs = crp_scs + 1;
                        dly = dly + crp_list(idx,1) - crp_list(idx,2) + 1;
                        blg = blg - 1;
                        crp_min_t = crp_list(idx,1) + 1;
                    else
                        crp_trans_idx = crp_list(:,3) == 1;
                        crp_trans_list = crp_list(crp_trans_idx,:);
                        crp_min_t = crp_trans_list(1,1) + 1;
                        crp_sect = sum(crp_trans_list(:,1) < crp_min_t);
                        if crp_sect == 1
                            scs = scs + 1;
                            crp_scs = crp_scs + 1;
                            dly = dly + crp_trans_list(1,1) - crp_trans_list(1,2) + 1;
                            crp_trans_list(1,3) = -1;
                            blg = blg - 1;
                            crp_trans_list(2:end,1) = crp_trans_list(1,1) + 1;
                        else
                            crp_sect_num = crp_sect - 1;
                            while crp_sect_num > 0 && crp_min_t < ENDTIME
                                crp_min_t = crp_trans_list(crp_sect,1) + 1;
                                crp_sect = sum(crp_trans_list(:,1) < crp_min_t);
                                crp_trans_list(1:crp_sect,3) = 2;
                                if crp_sect == crp_sect_num + 1
                                    break;
                                else
                                    crp_sect_num = crp_sect - 1;
                                end
                            end
                            crp_trans_list(1:crp_sect,1) = crp_min_t ...
                                + exprnd(1/crp_mu,crp_sect,1);
                            crp_trans_list(:,3) = crp_trans_list(:,3) - 1;  % particapate: 1, not: 0
                            crp_trans_list = sortrows(crp_trans_list,1);
                        end
                        crp_list(crp_trans_idx,:) = crp_trans_list;
                    end
                end
                crp_t = crp_t + crp_min_t - crp_start_t;
                
                if sum(crp_list(:,3) < -1) > 0
                    disp FALSE_CRP_SCS
                end
                % CRP end
                if sum(crp_list(:,3) >= 0) > 0
                    crp_list(crp_list(:,3)>=0,1) = crp_min_t + exprnd(1/mu,sum(crp_list(:,3)>=0),1);
                    crp_list(crp_list(:,3)>=0,3) = 1;
                    crp_list = sortrows(crp_list,1);
                end
                pkt_list(ptr:ptr+blg_end,:) = crp_list;
                channel_mint(flag) = crp_min_t;
            end
            
            ptr = scs + 1;
            pkt_list(ptr:scs+blg,:) = sortrows(pkt_list(ptr:scs+blg,:),1);
            
            % New arrival in channel
            if flag_status(channel_flag) == 0
                crp_min_t = channel_mint(channel_mint == min(channel_mint));
                while sum(pkt_list(ptr:end,1) < crp_min_t) > 0
                    % ac_blg = sum(pkt_list(:,2) < min_t) - scs;
                    % mu = 0.6468 / max(ac_blg,1);
                    cnt = cnt + 1;
                    if blg == 0
                        pkt_list(ptr,1) = pkt_list(ptr,1) + 0;
                        pkt_list(ptr,3) = 1;    % stack in backlog list
                        blg = 1;
                    end
                    min_t_temp = pkt_list(ptr,1) + 1;    % packet length equals 1
                    new_pkt = sum(pkt_list(ptr+blg:end,1) < min_t_temp);
                    if new_pkt > 0
                        bof = 0;
                        min_t_temp = min(min(pkt_list(ptr+blg:scs+blg+new_pkt,1)+bof)+1, min_t_temp);
                        new_blg = sum(pkt_list(ptr+blg:scs+blg+new_pkt,1) < min_t_temp);
                        pkt_list(ptr+blg:scs+blg+new_blg,1) = pkt_list(ptr+blg:scs+blg+new_blg,1) + bof;
                        pkt_list(ptr+blg:scs+blg+new_blg,3) = 1;
                        blg = blg + new_blg;
                        pkt_list(ptr:scs+blg,:) = sortrows(pkt_list(ptr:scs+blg,:),1);
                    end
                    idle_t = idle_t + pkt_list(ptr,1) - min_t;
                    if idle_t < 0
                        disp FALSE_IDLE_MINUS
                        % return
                    end
                    min_t = pkt_list(ptr,1) + 1;
                    sect = sum(pkt_list(ptr:ptr+blg-1,1) < min_t);
                    if sect == 1
                        scs = scs + 1;
                        dly = dly + pkt_list(ptr,1) - pkt_list(ptr,2) + 1;
                        pkt_list(ptr,3) = -1;   % column 3 == -1 => scs
                        ptr = ptr + 1;
                        blg = blg - 1;
                    else
                        cnt_coll = cnt_coll + 1;
                        blg_end = sect - 1;
                        coll_start_t = pkt_list(ptr,1);
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
                        coll_t = coll_t + min_t - coll_start_t;
                        pkt_list(ptr:ptr+blg_end,1) = min_t + exprnd(1/mu,sect,1);
                        pkt_list(ptr:scs+blg,:) = sortrows(pkt_list(ptr:scs+blg,:),1);
                    end
                    if sum(pkt_list(:,2) < min_t) - scs ~= blg
                        disp FALSE_BLG2
                        % return
                    end
                end
                channel_flag(channel_mint == min(channel_mint)) = 0;
            end
            pkt_list(ptr:scs+blg,:) = sortrows(pkt_list(ptr:scs+blg,:),1);
            if sum(pkt_list(:,2) < min_t) - scs ~= blg
                disp FALSE_BLG1
                % return
            end
        end
        ac_blg = sum(pkt_list(:,2) < min_t) - scs;
        if ac_blg ~= blg
            disp FALSE_BLG
            % return
        end
        prev_end_t = min_t;
        mu = 0.6468 / max(ac_blg,1); % a temp setting
    end
    thrpt_list(ldx) = scs / min_t;
    dly_list(ldx) = dly / scs;
    avg_coll(ldx) = coll_t / cnt_coll;
    avg_idle(ldx) = idle_t / cnt;
    % crp_scs_list(ldx) = crp_scs / crp_cnt;
    crp_len(ldx) = crp_t / crp_cnt;
    crp_avg_coll(ldx) = crp_coll / crp_cnt;
end
toc

figure
plot(lambda,thrpt_list,'LineWidth',1)
legend('Main Channel Throughput','Location','southeast','Interpreter','latex','FontSize',14.4)
grid on
% xlim([0 0.36])
xlabel('$\lambda$','Interpreter','latex','FontSize',17.6)
ylabel('Throughput (packet/sec)','Interpreter','latex','FontSize',17.6)
title('Pure ALOHA CRP','Interpreter','latex','FontSize',17.6)

figure
plot(lambda,avg_coll,lambda,avg_idle,lambda,avg_idle + avg_coll,lambda,crp_len,'LineWidth',1)
legend('Average Collision','Average Idle','Average Idle + Coll','Average CRP Length','Location','northeast','Interpreter','latex','FontSize',14.4)
grid on
% xlim([0 0.36])
xlabel('$\lambda$','Interpreter','latex','FontSize',17.6)
ylabel('Idle Time (sec)','Interpreter','latex','FontSize',17.6)
title('Pure ALOHA CRP','Interpreter','latex','FontSize',17.6)

figure
plot(lambda,crp_len ./ (avg_coll + avg_idle),'LineWidth',1)
legend('CRP:Trans','Location','southeast','Interpreter','latex','FontSize',14.4)
grid on
% xlim([0 0.36])
xlabel('$\lambda$','Interpreter','latex','FontSize',17.6)
ylabel('CRP:(Idle+Collision)','Interpreter','latex','FontSize',17.6)
title('Pure ALOHA CRP','Interpreter','latex','FontSize',17.6)