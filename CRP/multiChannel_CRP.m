clear

THEATA = 0.99;
ENDTIME = 2e4;
CHANNEL = 2;

lambda = 0.01:0.01:0.4;

thrpt_list = zeros(length(lambda),1);
dly_list = zeros(length(lambda),1);

crp_thrpt = zeros(length(lambda),1);

crp_mu = 0.5 / 2;
lambda_recur = 0.2;

tic
for ldx = 1:length(lambda)
    % initialize param & packet list
    num = ceil(2 * lambda(ldx) * ENDTIME);
    ptr = ones(CHANNEL,1);
    blg = zeros(CHANNEL,1);
    scs = blg;
    dly = blg;
    cnt = blg;
    min_t = blg;
    status = blg;
    blg_end = blg;
    mu = blg;
    mu(1:2) = 0.6468 / 10;

    crp_flag = 0;
    wait_flag = blg;

    crp_scs = 0;

    pkt_list = zeros(num,3,CHANNEL);
    for idx = 1:CHANNEL
        pkt_list(:,1,idx) = cumsum(exprnd(1/lambda(ldx), num, 1));
        pkt_list(:,2,idx) = pkt_list(:,1,idx);
    end
    
    while min_t(1) < ENDTIME || min_t(2) < ENDTIME
        
        for idx = 1:2
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
            % idle_t = idle_t + pkt_list(ptr(idx),1) - min_t;
            % if idle_t < 0
            %     disp FALSE_IDLE_MINUS
            %     % return
            % end
            min_t(idx) = pkt_list(ptr(idx),1,idx) + 1;
            sect = sum(pkt_list(ptr(idx):scs(idx)+blg(idx),1,idx) < min_t(idx));
            if sect == 1
                scs(idx) = scs(idx) + 1;
                dly(idx) = dly(idx) + pkt_list(ptr(idx),1,idx) - pkt_list(ptr(idx),2,idx) + 1;
                pkt_list(ptr(idx),3,idx) = -1;   % column 3 == -1 => scs(idx)
                ptr(idx) = ptr(idx) + 1;
                blg(idx) = blg(idx) - 1;
                status(idx) = 1;
            else
                coll_start_t = pkt_list(ptr(idx),1,idx);
                blg_end(idx) = sect - 1;
                % collision period
                while blg_end(idx) > 0 && min_t(idx) < ENDTIME
                    min_t(idx) = pkt_list(ptr(idx)+blg_end(idx),1,idx) + 1;
                    new_blg = sum(pkt_list(ptr(idx)+blg(idx):end,1,idx) < min_t(idx));
                    if new_blg > 0
                        pkt_list(ptr(idx)+blg(idx):scs(idx)+blg(idx)+new_blg,1,idx) = pkt_list( ...
                            ptr(idx)+blg(idx):scs(idx)+blg(idx)+new_blg,1,idx) + exprnd(1/mu(idx),new_blg,1);
                        pkt_list(ptr(idx)+blg(idx):scs(idx)+blg(idx)+new_blg,3,idx) = 1;
                        blg(idx) = blg(idx) + new_blg;
                        pkt_list(ptr(idx):scs(idx)+blg(idx),:,idx) = sortrows(pkt_list(ptr(idx):scs(idx)+blg(idx),:,idx),1);
                    end
                    sect = sum(pkt_list(ptr(idx):scs(idx)+blg(idx),1,idx) < min_t(idx));
                    if sect == blg_end(idx) + 1
                        break;
                    else
                        blg_end(idx) = sect - 1;
                    end
                end
                status(idx) = 2;
            end
        end
        % CRP procedure
        if sum(status == 2) == 1
            crp_idx = find(status == 2);
            crp_flag = 1;
        elseif sum(status == 2) == 2
            if min_t(1) <= min_t(2)
                crp_idx = 1;
            else
                crp_idx = 2;
            end
            crp_flag = 1;
            jdx = 3 - crp_idx;
            pkt_list(ptr(jdx):ptr(jdx)+blg_end(jdx),1,jdx) = min_t(jdx) + exprnd(1/mu(jdx),blg_end(jdx)+1,1);
            pkt_list(ptr(jdx):scs(jdx)+blg(jdx),:,jdx) = sortrows(pkt_list(ptr(jdx):scs(jdx)+blg(jdx),:,jdx),1);
        end
        if crp_flag == 1
            pkt_list(ptr(crp_idx):ptr(crp_idx)+blg_end(crp_idx),1,crp_idx) = min_t(crp_idx);
            crp_start_t = min_t(crp_idx);
            crp_list = pkt_list(ptr(crp_idx):ptr(crp_idx)+blg_end(crp_idx),:,crp_idx);
            if sum(crp_list(:,3) == 1) <= 1
                disp FALSE_COLL
            end
            crp_min_t = 0;
            crp_blg = sum(crp_list(:,3)==1);
            crp_list(:,1) = crp_list(:,1) + exprnd(1/crp_mu,crp_blg,1);
            crp_list = sortrows(crp_list,1);
            while sum(crp_list(:,3) == 1)  > 0 && crp_min_t < ENDTIME
                if sum(crp_list(:,3) == 1) == 1
                    idx = find(crp_list(:,3) == 1);
                    crp_list(idx,3) = -1;   % column 3 == -1 => scs
                    scs(crp_idx) = scs(crp_idx) + 1;
                    crp_scs = crp_scs + 1;
                    dly(crp_idx) = dly(crp_idx) + crp_list(idx,1) - crp_list(idx,2) + 1;
                    blg(crp_idx) = blg(crp_idx) - 1;
                    crp_min_t = crp_list(idx,1) + 1;
                else
                    crp_trans_idx = crp_list(:,3) == 1;
                    crp_trans_list = crp_list(crp_trans_idx,:);
                    crp_min_t = crp_trans_list(1,1) + 1;
                    crp_sect = sum(crp_trans_list(:,1) < crp_min_t);
                    if crp_sect == 1
                        scs(crp_idx) = scs(crp_idx) + 1;
                        crp_scs = crp_scs + 1;
                        dly(crp_idx) = dly(crp_idx) + crp_trans_list(1,1) - crp_trans_list(1,2) + 1;
                        crp_trans_list(1,3) = -1;
                        blg(crp_idx) = blg(crp_idx) - 1;
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

            if sum(crp_list(:,3) < -1) > 0
                disp FALSE_CRP_SCS
            end

            % CRP end
            if sum(crp_list(:,3) >= 0) > 0
                crp_list(crp_list(:,3)>=0,1) = crp_min_t + exprnd(1/mu(crp_idx),sum(crp_list(:,3)>=0),1);
                crp_list(crp_list(:,3)>=0,3) = 1;
                crp_list = sortrows(crp_list,1);
            end
            pkt_list(ptr(crp_idx):ptr(crp_idx)+blg_end(crp_idx),:,crp_idx) = crp_list;
            ptr(crp_idx) = scs(crp_idx) + 1;
            pkt_list(ptr(crp_idx):scs(crp_idx)+blg(crp_idx),:,crp_idx) = sortrows(pkt_list(ptr(crp_idx):scs(crp_idx)+blg(crp_idx),:,crp_idx),1);

            %% 
            for idx = 1:2
                while sum(pkt_list(ptr(idx):end,1,idx) < crp_min_t) > 0 && min_t(idx) < ENDTIME
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
                    % idle_t = idle_t + pkt_list(ptr(idx),1) - min_t;
                    % if idle_t < 0
                    %     disp FALSE_IDLE_MINUS
                    %     % return
                    % end
                    min_t(idx) = pkt_list(ptr(idx),1,idx) + 1;
                    sect = sum(pkt_list(ptr(idx):scs(idx)+blg(idx),1,idx) < min_t(idx));
                    if sect == 1
                        scs(idx) = scs(idx) + 1;
                        dly(idx) = dly(idx) + pkt_list(ptr(idx),1,idx) - pkt_list(ptr(idx),2,idx) + 1;
                        pkt_list(ptr(idx),3,idx) = -1;   % column 3 == -1 => scs(idx)
                        ptr(idx) = ptr(idx) + 1;
                        blg(idx) = blg(idx) - 1;
                    else
                        coll_start_t = pkt_list(ptr(idx),1,idx);
                        blg_end(idx) = sect - 1;
                        % collision period
                        while blg_end(idx) > 0 && min_t(idx) < ENDTIME
                            min_t(idx) = pkt_list(ptr(idx)+blg_end(idx),1,idx) + 1;
                            new_blg = sum(pkt_list(ptr(idx)+blg(idx):end,1,idx) < min_t(idx));
                            if new_blg > 0
                                pkt_list(ptr(idx)+blg(idx):scs(idx)+blg(idx)+new_blg,1,idx) = pkt_list( ...
                                    ptr(idx)+blg(idx):scs(idx)+blg(idx)+new_blg,1,idx) + exprnd(1/mu(idx),new_blg,1);
                                pkt_list(ptr(idx)+blg(idx):scs(idx)+blg(idx)+new_blg,3,idx) = 1;
                                blg(idx) = blg(idx) + new_blg;
                                pkt_list(ptr(idx):scs(idx)+blg(idx),:,idx) = sortrows(pkt_list(ptr(idx):scs(idx)+blg(idx),:,idx),1);
                            end
                            sect = sum(pkt_list(ptr(idx):scs(idx)+blg(idx),1,idx) < min_t(idx));
                            if sect == blg_end(idx) + 1
                                break;
                            else
                                blg_end(idx) = sect - 1;
                            end
                        end
                        pkt_list(ptr(idx):ptr(idx)+blg_end(idx),1,idx) = min_t(idx) + exprnd(1/mu(idx),blg_end(idx)+1,1);
                        pkt_list(ptr(idx):scs(idx)+blg(idx),:,idx) = sortrows(pkt_list(ptr(idx):scs(idx)+blg(idx),:,idx),1);
                    end
                end
            end
            crp_flag = 0;
        end
        mu = 0.6468 ./ max(blg,1);
    end
    thrpt_list(ldx) = sum(scs ./ min_t) / 3;
    dly_list(ldx) = sum(dly ./ scs) / 2;
    crp_thrpt(ldx) = crp_scs / sum(min_t) * 2;
end
toc

figure
plot(lambda,crp_thrpt,lambda,thrpt_list,'LineWidth',1)
legend('CRP Thrpughput','Total Throughput','Location','southeast','Interpreter','latex','FontSize',14.4)
grid on
% xlim([0 0.36])
xlabel('$\lambda$','Interpreter','latex','FontSize',17.6)
ylabel('Throughput (packet/sec)','Interpreter','latex','FontSize',17.6)
title('Pure ALOHA CRP','Interpreter','latex','FontSize',17.6)

figure
plot(lambda,dly_list,'LineWidth',1)
legend('Delay','Location','southeast','Interpreter','latex','FontSize',14.4)
grid on
% xlim([0 0.36])
xlabel('$\lambda$','Interpreter','latex','FontSize',17.6)
ylabel('Delay (sec)','Interpreter','latex','FontSize',17.6)
title('Pure ALOHA CRP','Interpreter','latex','FontSize',17.6)