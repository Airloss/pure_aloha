clear

THEATA = 0.99;
ENDTIME = 2e4;
CHANNEL = 2;

lambda = 0.02:0.02:0.6;
betaT = 0.6;

thrpt_list = zeros(length(lambda),1);
chanl_thrpt = zeros(length(lambda),1);
dly_list = zeros(length(lambda),1);

crp_thrpt = zeros(length(lambda),1);
crp_thrpt_ideal = zeros(length(lambda),1);
crp_len = zeros(length(lambda),1);
crp_invov = zeros(length(lambda),1);

tic
parfor (ldx = 1:length(lambda),6)
    % initialize param & packet list
    num = ceil(1.5 * lambda(ldx) * ENDTIME);
    ptr = ones(CHANNEL,1);
    blg = zeros(CHANNEL,1);
    scs = blg;
    dly = blg;
    cnt = blg;
    min_t = blg;
    status = blg;
    blg_end = blg;
    mu = blg;
    mu(1:2) = 0.6468 / 2;

    crp_flag = 0;
    wait_flag = blg;

    crp_t = 0;
    crp_cnt = 0;
    crp_scs = 0;
    crp_idx = 0;
    crp_num = 0;

    pkt_list = zeros(num,3,CHANNEL);
    for i1 = 1:CHANNEL
        pkt_list(:,1,i1) = cumsum(exprnd(1/lambda(ldx), num, 1));
        pkt_list(:,2,i1) = pkt_list(:,1,i1);
    end
    
    while min_t(1) < ENDTIME || min_t(2) < ENDTIME
        %% Before Transmission
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
                        pkt_list(ptr(idx):scs(idx)+blg(idx),:,idx) = ...
                            sortrows(pkt_list(ptr(idx):scs(idx)+blg(idx),:,idx),1);
                    end
                    sect = sum(pkt_list(ptr(idx):scs(idx)+blg(idx),1,idx) < min_t(idx));
                    if sect == blg_end(idx) + 1
                        break;
                    else
%                         disp(1);
%                         disp(pkt_list(ptr(idx):ptr(idx)+blg_end(idx),1,idx)')
                        blg_end(idx) = sect - 1;
                    end
                end
                status(idx) = 2;
            end
        end
        %% CRP procedure
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
            crp_cnt = crp_cnt + 1;
            crp_scs = crp_scs + 2;
            crp_min_t = min_t(crp_idx) + 2;
            crp_t = crp_t + 2;
            crp_list = pkt_list(ptr(crp_idx):ptr(crp_idx)+blg_end(crp_idx),:,crp_idx);
            % first pkt success
            scs(crp_idx) = scs(crp_idx) + 1;
            dly(crp_idx) = dly(crp_idx) + pkt_list(ptr(crp_idx),1,crp_idx) - pkt_list(ptr(crp_idx),2,crp_idx) + 1;
            pkt_list(ptr(crp_idx),3,crp_idx) = -1;
            crp_list(1,3) = -1;
            % last pkt success
            scs(crp_idx) = scs(crp_idx) + 1;
            dly(crp_idx) = dly(crp_idx) + pkt_list(ptr(crp_idx)+blg_end(crp_idx),1,crp_idx) - pkt_list(ptr(crp_idx)+blg_end(crp_idx),2,crp_idx) + 1;
            pkt_list(ptr(crp_idx)+blg_end(crp_idx),3,crp_idx) = -1;
            crp_list(end,3) = -1;
            crp_list(2:end-1,1) = min_t(crp_idx) + 2 + exprnd(1/mu(crp_idx),blg_end(crp_idx)-1,1);
            crp_list = sortrows(crp_list,1);
            % end
            ptr(crp_idx) = ptr(crp_idx) + 2;
            blg(crp_idx) = blg(crp_idx) - 2;
            pkt_list(ptr(crp_idx)+1:ptr(crp_idx)+blg_end(crp_idx)-1,1,crp_idx) = 2 + min_t(crp_idx) + exprnd(1/mu(crp_idx),blg_end(crp_idx)-1,1);
            pkt_list(ptr(crp_idx):ptr(crp_idx)+blg(crp_idx),:,crp_idx) = sortrows(pkt_list(ptr(crp_idx):ptr(crp_idx)+blg(crp_idx),:,crp_idx),1);

            
            %% During CRP
            for ii = 1:2
                while sum(pkt_list(ptr(ii):end,1,ii) < crp_min_t) > 0
                    if blg(ii) == 0
                        pkt_list(ptr(ii),1,ii) = pkt_list(ptr(ii),1,ii) + exprnd(mu(ii),1);
                        pkt_list(ptr(ii),3,ii) = 1;    % stack in backlog list
                        blg(ii) = 1;
                    end
                    min_t_temp = pkt_list(ptr(ii),1,ii) + 1;    % packet length equals 1
                    new_pkt = sum(pkt_list(ptr(ii)+blg(ii):end,1,ii) < min_t_temp);
                    if new_pkt > 0
                        bof = exprnd(1/mu(ii), new_pkt, 1);
                        min_t_temp = min(min( ...
                            pkt_list(ptr(ii)+blg(ii):scs(ii)+blg(ii)+new_pkt,1,ii)+bof)+1, min_t_temp);
                        new_blg = sum(pkt_list(ptr(ii)+blg(ii):scs(ii)+blg(ii)+new_pkt,1,ii) < min_t_temp);
                        pkt_list(ptr(ii)+blg(ii):scs(ii)+blg(ii)+new_blg,1,ii) = pkt_list( ...
                            ptr(ii)+blg(ii):scs(ii)+blg(ii)+new_blg,1,ii) + bof(1:new_blg);
                        pkt_list(ptr(ii)+blg(ii):scs(ii)+blg(ii)+new_blg,3,ii) = 1;
                        blg(ii) = blg(ii) + new_blg;
                        pkt_list(ptr(ii):scs(ii)+blg(ii),:,ii) = sortrows(pkt_list(ptr(ii):scs(ii)+blg(ii),:,ii),1);
                    end
                    % idle_t = idle_t + pkt_list(ptr(ii),1) - min_t;
                    % if idle_t < 0
                    %     disp FALSE_IDLE_MINUS
                    %     % return
                    % end
                    min_t(ii) = pkt_list(ptr(ii),1,ii) + 1;
                    sect = sum(pkt_list(ptr(ii):scs(ii)+blg(ii),1,ii) < min_t(ii));
                    if sect == 1
                        scs(ii) = scs(ii) + 1;
                        dly(ii) = dly(ii) + pkt_list(ptr(ii),1,ii) - pkt_list(ptr(ii),2,ii) + 1;
                        pkt_list(ptr(ii),3,ii) = -1;   % column 3 == -1 => scs(ii)
                        ptr(ii) = ptr(ii) + 1;
                        blg(ii) = blg(ii) - 1;
                        status(ii) = 1;
                    else
                        coll_start_t = pkt_list(ptr(ii),1,ii);
                        blg_end(ii) = sect - 1;
                        % collision period
                        while blg_end(ii) > 0 && min_t(ii) < ENDTIME
                            min_t(ii) = pkt_list(ptr(ii)+blg_end(ii),1,ii) + 1;
                            new_blg = sum(pkt_list(ptr(ii)+blg(ii):end,1,ii) < min_t(ii));
                            if new_blg > 0
                                pkt_list(ptr(ii)+blg(ii):scs(ii)+blg(ii)+new_blg,1,ii) = pkt_list( ...
                                    ptr(ii)+blg(ii):scs(ii)+blg(ii)+new_blg,1,ii) + exprnd(1/mu(ii),new_blg,1);
                                pkt_list(ptr(ii)+blg(ii):scs(ii)+blg(ii)+new_blg,3,ii) = 1;
                                blg(ii) = blg(ii) + new_blg;
                                pkt_list(ptr(ii):scs(ii)+blg(ii),:,ii) = sortrows(pkt_list(ptr(ii):scs(ii)+blg(ii),:,ii),1);
                            end
                            sect = sum(pkt_list(ptr(ii):scs(ii)+blg(ii),1,ii) < min_t(ii));
                            if sect == blg_end(ii) + 1
                                break;
                            else
%                                 disp(2);
%                                 disp(pkt_list(ptr(idx):ptr(idx)+blg_end(idx),1,idx)')
                                blg_end(ii) = sect - 1;
                            end
                        end
                        status(ii) = 2;
                        if min_t(ii) > crp_min_t
                            break;
                        else
                            pkt_list(ptr(ii):ptr(ii)+blg_end(ii),1,ii) = min_t(ii) + exprnd(1/mu(ii),blg_end(ii)+1,1);
                            pkt_list(ptr(ii):scs(ii)+blg(ii),:,ii) = sortrows(pkt_list(ptr(ii):scs(ii)+blg(ii),:,ii),1);
                        end
                    end
                end
            end
            crp_flag = 0;
        end
        mu = betaT ./ max(blg,1);
    end
    thrpt_list(ldx) = sum(scs ./ min_t) / 3;
    chanl_thrpt(ldx) = (sum(scs)-crp_scs) / sum(min_t);
    dly_list(ldx) = sum(dly ./ scs) / 2;
    crp_thrpt(ldx) = crp_scs / sum(min_t) * 2;
    crp_thrpt_ideal(ldx) = crp_scs / crp_t;
    crp_len(ldx) = crp_t / sum(min_t) * 2;
    % crp_invov(ldx) = crp_num / crp_cnt;
end
toc

pt = find(thrpt_list == max(thrpt_list));
disp(thrpt_list(pt));
disp(lambda(pt));

yy = ones(length(lambda),1);
yy = yy .* 0.184;

figure
plot(lambda,crp_thrpt,lambda,chanl_thrpt,lambda,thrpt_list,lambda,yy,'--','LineWidth',1.5)
legend('CRP Thrpughput','Channel Throughput','Total Throughput','Location','southeast','Interpreter','latex','FontSize',14.4)
grid on
% xlim([0 0.36])
xlabel('$\lambda$','Interpreter','latex','FontSize',17.6)
ylabel('Throughput (packet/sec)','Interpreter','latex','FontSize',17.6)
title('Pure ALOHA CRP','Interpreter','latex','FontSize',17.6)

figure
plot(lambda,dly_list,'LineWidth',1)
legend('Delay','Location','northwest','Interpreter','latex','FontSize',14.4)
grid on
ylim([0 300])
xlabel('$\lambda$','Interpreter','latex','FontSize',17.6)
ylabel('Delay (sec)','Interpreter','latex','FontSize',17.6)
title('Pure ALOHA CRP','Interpreter','latex','FontSize',17.6)

figure
plot(lambda,crp_len,'LineWidth',1.5)
legend('CRP Propotion','Location','southeast','Interpreter','latex','FontSize',14.4)
grid on
% xlim([0 0.36])
xlabel('$\lambda$','Interpreter','latex','FontSize',17.6)
ylabel('Channel Efficiency','Interpreter','latex','FontSize',17.6)
title('Pure ALOHA CRP','Interpreter','latex','FontSize',17.6)

figure
plot(lambda,crp_thrpt_ideal,'LineWidth',1.5)
legend('CRP Throughput Ideal','Location','southeast','Interpreter','latex','FontSize',14.4)
grid on
% xlim([0 0.36])
xlabel('$\lambda$','Interpreter','latex','FontSize',17.6)
ylabel('Channel Efficiency','Interpreter','latex','FontSize',17.6)
title('Pure ALOHA CRP','Interpreter','latex','FontSize',17.6)