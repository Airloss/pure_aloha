clear

THEATA = 0.99;
ENDTIME = 2e5;
CHANNEL = 2;

lambda = 0.03:0.01:0.36;
betaT = 0.5;
crp_mu = 0.25;

sys_thrpt = zeros(length(lambda),1);
chanl_thrpt = zeros(length(lambda),1);
dly_list = zeros(length(lambda),1);
blg_list = zeros(length(lambda),1);

crp_thrpt = zeros(length(lambda),1);
crp_thrpt_ideal = zeros(length(lambda),1);
crp_chnl_util = zeros(length(lambda),1);
crp_invov = zeros(length(lambda),1);

tic
parfor (ldx = 1:length(lambda),8)
    % initialize param & packet list
    num = ceil(1.5 * lambda(ldx) * ENDTIME);
    ptr = ones(CHANNEL,1);
    blg = zeros(CHANNEL,1);
    scs = blg;
    dly = blg;
    cnt = blg;
    min_t = blg;
    idle_t = blg;
    prev_end_t = blg;
    status = blg;
    blg_end = blg;
    l_recur = ptr .* 0.2;
    mu = blg;
    es_blg = blg;
    mu(1:end) = betaT / 2;
    es_blg(1:end) = 1;
    chnl_scs = scs;

    blg_diff = blg;
    crp_flag = 0;

    crp_t = 0;
    crp_cnt = 0;
    crp_scs = 0;
    crp_idx = 0;
    crp_num = 0;

    pkt_list = zeros(num,3,CHANNEL);
    for ii = 1:CHANNEL
        pkt_list(:,1,ii) = cumsum(exprnd(1/lambda(ldx), num, 1));
        pkt_list(:,2,ii) = pkt_list(:,1,ii);
    end
    
    while sum(min_t <= ENDTIME) > 0
        %% Before Transmission
        for idx = 1:CHANNEL
            cnt(idx) = cnt(idx) + 1;
            if blg(idx) == 0
                pkt_list(ptr(idx),1,idx) = pkt_list(ptr(idx),1,idx) + exprnd(mu(idx),1);
                pkt_list(ptr(idx),3,idx) = 1;    % stack in backlog list
                blg(idx) = 1;
            end
            min_t_temp = pkt_list(ptr(idx),1,idx) + 1;    % packet length equals 1
            new_pkt = sum(pkt_list(ptr(idx)+blg(idx):end,1,idx) < min_t_temp);
            if new_pkt > 0
                bof = exprnd(1/mu(idx), new_pkt, 1);
                min_t_temp = min(min( ...
                    pkt_list(ptr(idx)+blg(idx):scs(idx)+blg(idx)+new_pkt,1,idx)+bof)+1, min_t_temp);
                new_blg = sum(pkt_list(ptr(idx)+blg(idx):scs(idx)+blg(idx)+new_pkt,1,idx) < min_t_temp);
                pkt_list(ptr(idx)+blg(idx):scs(idx)+blg(idx)+new_blg,1,idx) = pkt_list( ...
                    ptr(idx)+blg(idx):scs(idx)+blg(idx)+new_blg,1,idx) + bof(1:new_blg);
                pkt_list(ptr(idx)+blg(idx):scs(idx)+blg(idx)+new_blg,3,idx) = 1;
                blg(idx) = blg(idx) + new_blg;
                pkt_list(ptr(idx):scs(idx)+blg(idx),:,idx) = sortrows(pkt_list(ptr(idx):scs(idx)+blg(idx),:,idx),1);
            end
            idle_t(idx) = pkt_list(ptr(idx),1,idx) - min_t(idx);
            if idle_t(idx) < 0
                disp FALSE_IDLE_T_1
            end
            min_t(idx) = pkt_list(ptr(idx),1,idx) + 1;
            sect = sum(pkt_list(ptr(idx):scs(idx)+blg(idx),1,idx) < min_t(idx));
            if sect == 1
                l_recur(idx) = THEATA * l_recur(idx) + (1 - THEATA) / (min_t(idx) - prev_end_t(idx));
                if l_recur(idx) < 0
                    disp FALSE_L_RECUR_1
                end
                es_blg(idx) = max(es_blg(idx) * exp(-mu(idx) * idle_t(idx)) + l_recur(idx),1);
                blg_diff(idx) = blg_diff(idx) + (es_blg(idx) - blg(idx));
                scs(idx) = scs(idx) + 1;
                chnl_scs(idx) = chnl_scs(idx) + 1;
                dly(idx) = dly(idx) + pkt_list(ptr(idx),1,idx) - pkt_list(ptr(idx),2,idx) + 1;
                pkt_list(ptr(idx),3,idx) = -1;   % column 3 == -1 => scs(idx)
                ptr(idx) = ptr(idx) + 1;
                blg(idx) = blg(idx) - 1;
                status(idx) = -1;
            else
                coll_start_t = pkt_list(ptr(idx),1,idx);
                blg_end(idx) = sect - 1;
                if sect < 1
                    disp FALSE_SECT_1
                end
                % collision period
                while blg_end(idx) > 0 && min_t(idx) < ENDTIME
                    min_t(idx) = pkt_list(ptr(idx)+blg_end(idx),1,idx) + 1;
                    new_blg = sum(pkt_list(ptr(idx)+blg(idx):end,1,idx) < min_t(idx));
                    if new_blg > 0
                        pkt_list(ptr(idx)+blg(idx):scs(idx)+blg(idx)+new_blg,1,idx) = pkt_list( ...
                            ptr(idx)+blg(idx):scs(idx)+blg(idx)+new_blg,1,idx) + exprnd(1/mu(idx),new_blg,1);
                        pkt_list(ptr(idx)+blg(idx):scs(idx)+blg(idx)+new_blg,3,idx) = 1;
                        blg(idx) = blg(idx) + new_blg;
                        pkt_list(ptr(idx):scs(idx)+blg(idx),:,idx) = ...
                            sortrows(pkt_list(ptr(idx):scs(idx)+blg(idx),:,idx),1);
                    end
                    sect = sum(pkt_list(ptr(idx):scs(idx)+blg(idx),1,idx) < min_t(idx));
                    if sect == blg_end(idx) + 1
                        break;
                    else
                        blg_end(idx) = sect - 1;
                    end
                end
                if blg_end(idx) < 1
                    disp FLASE_BLGEND
                end
                status(idx) = 1;
                l_recur(idx) = THEATA * l_recur(idx);
                if l_recur(idx) < 0
                    disp FALSE_L_RECUR_2
                end
                es_blg(idx) = max(1 + es_blg(idx) * exp(-mu(idx) * idle_t(idx)) + l_recur(idx) * (min_t(idx) - coll_start_t),1);
                blg_diff(idx) = blg_diff(idx) + (es_blg(idx) - blg(idx));
            end
        end
        mu = betaT ./ es_blg;
        prev_end_t = min_t;

        % Transmission end
        status_ = status == 1;
        if sum(status_) == 1
            crp_idx = find(status == 1);
            crp_flag = 1;
        elseif sum(status_) > 1
            pkt_list_ = zeros(CHANNEL,1);
            for jj = 1:CHANNEL
                pkt_list_(jj) = pkt_list(ptr(jj),1,jj);
            end
            pkt_list_(~status_) = max(pkt_list_) + 1;
            crp_idx = find(pkt_list_ == min(pkt_list_));
            status_(crp_idx) = 0;
            for kk = 1:CHANNEL
                if status_(kk) > 0
                    pkt_list(ptr(kk):ptr(kk)+blg_end(kk),1,kk) = min_t(kk) + exprnd(1/mu(kk),blg_end(kk)+1,1);
                    pkt_list(ptr(kk):scs(kk)+blg(kk),:,kk) = sortrows(pkt_list(ptr(kk):scs(kk)+blg(kk),:,kk),1);
                end
            end
            crp_flag = 1;
        end

        %% CRP
        if crp_flag == 1
            crp_cnt = crp_cnt + 1;
            pkt_list(ptr(crp_idx):ptr(crp_idx)+blg_end(crp_idx),1,crp_idx) = min_t(crp_idx);
            crp_start_t = min_t(crp_idx);
            crp_list = pkt_list(ptr(crp_idx):ptr(crp_idx)+blg_end(crp_idx),:,crp_idx);
            if sum(crp_list(:,3) == 1) <= 1
                disp FALSE_COLL
            end
            crp_scs_ = 0;
            crp_min_t = 0;
            crp_blg = sum(crp_list(:,3) == 1);
            crp_num = crp_num + crp_blg;
            crp_list(:,1) = crp_list(:,1) + exprnd(1/crp_mu,crp_blg,1);
            crp_list = sortrows(crp_list,1);
            while sum(crp_list(:,3) == 1)  > 0 && crp_min_t < ENDTIME
                if sum(crp_list(:,3) == 1) == 1
                    ic = find(crp_list(:,3) == 1);
                    crp_list(ic,3) = -1;   % column 3 == -1 => scs
                    scs(crp_idx) = scs(crp_idx) + 1;
                    crp_scs = crp_scs + 1;
                    crp_scs_ = crp_scs_ + 1;
                    dly(crp_idx) = dly(crp_idx) + crp_list(ic,1) - crp_list(ic,2) + 1;
                    crp_min_t = crp_list(ic,1) + 1;
                else
                    crp_trans_idx = crp_list(:,3) == 1;
                    crp_trans_list = crp_list(crp_trans_idx,:);
                    crp_min_t = crp_trans_list(1,1) + 1;
                    crp_sect = sum(crp_trans_list(:,1) < crp_min_t);
                    if crp_sect == 1
                        scs(crp_idx) = scs(crp_idx) + 1;
                        crp_scs = crp_scs + 1;
                        crp_scs_ = crp_scs_ + 1;
                        dly(crp_idx) = dly(crp_idx) + crp_trans_list(1,1) - crp_trans_list(1,2) + 1;
                        crp_trans_list(1,3) = -1;
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
                        crp_trans_list(1:crp_sect,1) = crp_min_t + exprnd(1/crp_mu,crp_sect,1);
                        crp_trans_list(:,3) = crp_trans_list(:,3) - 1;  % particapate: 1, not: 0
                        crp_trans_list = sortrows(crp_trans_list,1);
                    end
                    crp_list(crp_trans_idx,:) = crp_trans_list;
                end
            end
            crp_t = crp_t + crp_min_t - crp_start_t;

            if sum(crp_list(:,3) < -1) > 0 || crp_scs_ ~= sum(crp_list(:,3) < 0)
                disp FALSE_CRP_SCS
            end

            % CRP end
            crp_blg_end = sum(crp_list(:,3) >= 0);
            blg(crp_idx) = blg(crp_idx) - sum(crp_list(:,3) < 0);
            if crp_blg_end > 0
                crp_list(crp_list(:,3)>=0,1) = crp_min_t;
                crp_list(crp_list(:,3)>=0,3) = 0;
                crp_list = sortrows(crp_list,1);
            end
            pkt_list(ptr(crp_idx):ptr(crp_idx)+blg_end(crp_idx),:,crp_idx) = crp_list;
            ptr(crp_idx) = scs(crp_idx) + 1;
            if crp_blg_end > 0
                if blg(crp_idx) - crp_blg_end > 0
                    pkt_list(ptr(crp_idx):scs(crp_idx)+blg(crp_idx)-crp_blg_end,:,crp_idx) = pkt_list(ptr(crp_idx)+crp_blg_end:scs(crp_idx)+blg(crp_idx),:,crp_idx);
                    pkt_list(ptr(crp_idx)+blg(crp_idx)-crp_blg_end:scs(crp_idx)+blg(crp_idx),:,crp_idx) = crp_list(crp_blg-crp_blg_end+1:crp_blg,:);
                    blg(crp_idx) = blg(crp_idx) - crp_blg_end;
                    pkt_list(ptr(crp_idx):scs(crp_idx)+blg(crp_idx),:,crp_idx) = sortrows(pkt_list(ptr(crp_idx):scs(crp_idx)+blg(crp_idx),:,crp_idx),1);
                    pkt_list(ptr(crp_idx)+blg(crp_idx):end,:,crp_idx) = sortrows(pkt_list(ptr(crp_idx)+blg(crp_idx):end,:,crp_idx),1);
                else
                    if blg(crp_idx) - crp_blg_end ~= 0
                        disp FALSE_CRP_BLG_END
                    end
                    blg(crp_idx) = blg(crp_idx) - crp_blg_end;
                    pkt_list(ptr(crp_idx)+blg(crp_idx):end,:,crp_idx) = sortrows(pkt_list(ptr(crp_idx)+blg(crp_idx):end,:,crp_idx),1);
                end
            else
                pkt_list(ptr(crp_idx):scs(crp_idx)+blg(crp_idx),:,crp_idx) = sortrows(pkt_list(ptr(crp_idx):scs(crp_idx)+blg(crp_idx),:,crp_idx),1);
            end

            %% During CRP
            for jdx = 1:CHANNEL
                prev_min_t = min_t(jdx);
                while sum(pkt_list(ptr(jdx):end,1,jdx) < crp_min_t) > 0
                    cnt(jdx) = cnt(jdx) + 1;
                    if blg(jdx) == 0
                        pkt_list(ptr(jdx),1,jdx) = pkt_list(ptr(jdx),1,jdx) + exprnd(mu(jdx),1);
                        pkt_list(ptr(jdx),3,jdx) = 1;    % stack in backlog list
                        blg(jdx) = 1;
                    end
                    min_t_temp = pkt_list(ptr(jdx),1,jdx) + 1;    % packet length equals 1
                    new_pkt = sum(pkt_list(ptr(jdx)+blg(jdx):end,1,jdx) < min_t_temp);
                    if new_pkt > 0
                        bof = exprnd(1/mu(jdx), new_pkt, 1);
                        min_t_temp = min(min( ...
                            pkt_list(ptr(jdx)+blg(jdx):scs(jdx)+blg(jdx)+new_pkt,1,jdx)+bof)+1, min_t_temp);
                        new_blg = sum(pkt_list(ptr(jdx)+blg(jdx):scs(jdx)+blg(jdx)+new_pkt,1,jdx) < min_t_temp);
                        pkt_list(ptr(jdx)+blg(jdx):scs(jdx)+blg(jdx)+new_blg,1,jdx) = pkt_list( ...
                            ptr(jdx)+blg(jdx):scs(jdx)+blg(jdx)+new_blg,1,jdx) + bof(1:new_blg);
                        pkt_list(ptr(jdx)+blg(jdx):scs(jdx)+blg(jdx)+new_blg,3,jdx) = 1;
                        blg(jdx) = blg(jdx) + new_blg;
                        pkt_list(ptr(jdx):scs(jdx)+blg(jdx),:,jdx) = sortrows(pkt_list(ptr(jdx):scs(jdx)+blg(jdx),:,jdx),1);
                    end
                    idle_t(jdx) = pkt_list(ptr(jdx),1,jdx) - min_t(jdx);
                    if idle_t(idx) < 0
                        disp FALSE_IDLE_T_2
                    end
                    min_t(jdx) = pkt_list(ptr(jdx),1,jdx) + 1;
                    if min_t(jdx) > crp_min_t
                        min_t(jdx) = prev_min_t;
                        break;
                    end
                    sect = sum(pkt_list(ptr(jdx):scs(jdx)+blg(jdx),1,jdx) < min_t(jdx));
                    if sect == 1
                        prev_min_t = min_t(jdx);
                        scs(jdx) = scs(jdx) + 1;
                        chnl_scs(jdx) = chnl_scs(jdx) + 1;
                        dly(jdx) = dly(jdx) + pkt_list(ptr(jdx),1,jdx) - pkt_list(ptr(jdx),2,jdx) + 1;
                        pkt_list(ptr(jdx),3,jdx) = -1;   % column 3 == -1 => scs(jdx)
                        ptr(jdx) = ptr(jdx) + 1;
                        blg(jdx) = blg(jdx) - 1;
                        status(jdx) = -1;
                        l_recur(jdx) = THEATA * l_recur(jdx) + (1 - THEATA) / (min_t(jdx) - prev_end_t(jdx));
                        if l_recur(jdx) < 0
                            disp FALSE_L_RECUR_3
                        end
                        es_blg(jdx) = max(es_blg(jdx) * exp(-mu(jdx) * idle_t(jdx)) + l_recur(jdx),1);
                        blg_diff(jdx) = blg_diff(jdx) + (es_blg(jdx) - blg(jdx));
                    else
                        coll_start_t = pkt_list(ptr(jdx),1,jdx);
                        blg_end(jdx) = sect - 1;
                        % collision period
                        while blg_end(jdx) > 0 && min_t(jdx) < ENDTIME
                            min_t(jdx) = pkt_list(ptr(jdx)+blg_end(jdx),1,jdx) + 1;
                            new_blg = sum(pkt_list(ptr(jdx)+blg(jdx):end,1,jdx) < min_t(jdx));
                            if new_blg > 0
                                pkt_list(ptr(jdx)+blg(jdx):scs(jdx)+blg(jdx)+new_blg,1,jdx) = pkt_list( ...
                                    ptr(jdx)+blg(jdx):scs(jdx)+blg(jdx)+new_blg,1,jdx) + exprnd(1/mu(jdx),new_blg,1);
                                pkt_list(ptr(jdx)+blg(jdx):scs(jdx)+blg(jdx)+new_blg,3,jdx) = 1;
                                blg(jdx) = blg(jdx) + new_blg;
                                pkt_list(ptr(jdx):scs(jdx)+blg(jdx),:,jdx) = sortrows(pkt_list(ptr(jdx):scs(jdx)+blg(jdx),:,jdx),1);
                            end
                            sect = sum(pkt_list(ptr(jdx):scs(jdx)+blg(jdx),1,jdx) < min_t(jdx));
                            if sect == blg_end(jdx) + 1
                                break;
                            else
                                blg_end(jdx) = sect - 1;
                            end
                        end
                        l_recur(jdx) = THEATA * l_recur(jdx);
                        if l_recur(jdx) < 0
                            disp FALSE_L_RECUR_4
                        end
                        es_blg(jdx) = max(1 + es_blg(jdx) * exp(-mu(jdx) * idle_t(jdx)) + l_recur(jdx) * (min_t(jdx) - coll_start_t),1);
                        blg_diff(jdx) = blg_diff(jdx) + (es_blg(jdx) - blg(jdx));
                        status(jdx) = 1;
                        if min_t(jdx) > crp_min_t
                            min_t(jdx) = prev_min_t;
                            break;
                        else
                            min_t(jdx) = prev_min_t;
                            pkt_list(ptr(jdx):ptr(jdx)+blg_end(jdx),1,jdx) = min_t(jdx) + exprnd(1/mu(jdx),blg_end(jdx)+1,1);
                            pkt_list(ptr(jdx):scs(jdx)+blg(jdx),:,jdx) = sortrows(pkt_list(ptr(jdx):scs(jdx)+blg(jdx),:,jdx),1);
                        end
                    end
                    mu(jdx) = betaT / es_blg(jdx);
                    prev_end_t(jdx) = min_t(jdx);
                end
            end
            crp_flag = 0;
        end
        prev_end_t = min_t;
        mu = betaT ./ es_blg;
    end
    if sum(chnl_scs) + crp_scs ~= sum(scs)
        disp(ldx);
    end
    sys_thrpt(ldx) = sum(scs ./ min_t) / (CHANNEL+1);
    chanl_thrpt(ldx) = mean(chnl_scs ./ min_t);
    dly_list(ldx) = mean(dly ./ scs);
    crp_thrpt(ldx) = crp_scs / mean(min_t);
    crp_thrpt_ideal(ldx) = crp_scs / crp_t;
    crp_chnl_util(ldx) = crp_t / mean(min_t);
    crp_invov(ldx) = crp_num / crp_cnt;
    blg_list(ldx) = mean(blg_diff ./ cnt);
end
toc

pt = find(sys_thrpt == max(sys_thrpt));
disp(sys_thrpt(pt));
disp(lambda(pt) * CHANNEL / (CHANNEL+1));

yy = ones(length(lambda),1);
yy = yy .* 0.184;

ftitle = sprintf('%d contention channels Unslotted ALOHA',CHANNEL);
xaxis_ = lambda * CHANNEL / (CHANNEL+1);

figure
plot(xaxis_,crp_thrpt,xaxis_,chanl_thrpt,xaxis_,sys_thrpt,xaxis_,yy,'--','LineWidth',1.5)
legend('CRP Thrpughput','Channel Throughput','Total Throughput','Location','southeast','Interpreter','latex','FontSize',14.4)
grid on
% xlim([0 0.36])
xlabel('$\lambda$','Interpreter','latex','FontSize',17.6)
ylabel('Throughput (packet/sec)','Interpreter','latex','FontSize',17.6)
title(ftitle,'Interpreter','latex','FontSize',17.6)

figure
plot(xaxis_,dly_list,'LineWidth',1)
legend('Delay','Location','northwest','Interpreter','latex','FontSize',14.4)
grid on
ylim([0 300])
xlabel('$\lambda$','Interpreter','latex','FontSize',17.6)
ylabel('Delay (sec)','Interpreter','latex','FontSize',17.6)
title('Unslotted ALOHA CRP','Interpreter','latex','FontSize',17.6)

figure
plot(xaxis_,crp_chnl_util,'LineWidth',1.5)
legend(ftitle,'Location','southeast','Interpreter','latex','FontSize',14.4)
grid on
% xlim([0 0.36])
xlabel('$\lambda$','Interpreter','latex','FontSize',17.6)
ylabel('Channel Utilization','Interpreter','latex','FontSize',17.6)
title('Unslotted ALOHA CRP','Interpreter','latex','FontSize',17.6)

figure
plot(xaxis_,blg_list,'LineWidth',1.5)
legend(ftitle,'Location','southeast','Interpreter','latex','FontSize',14.4)
grid on
% ylim([-100 10])
xlabel('$\lambda$','Interpreter','latex','FontSize',17.6)
ylabel('Backoff Difference','Interpreter','latex','FontSize',17.6)
title('Unslotted ALOHA CRP','Interpreter','latex','FontSize',17.6)