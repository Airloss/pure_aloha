clear

ENDTIME = 1e6;

lambda = 0.01:0.01:0.56;

sys_thrpt = zeros(length(lambda),1);
chanl_thrpt = zeros(length(lambda),1);
dly_list = zeros(length(lambda),1);

crp_thrpt = zeros(length(lambda),1);
crp_thrpt_ideal = zeros(length(lambda),1);
crp_chnl_util = zeros(length(lambda),1);
crp_invov = zeros(length(lambda),1);

tic
parfor ldx = 1:length(lambda)
    num = ceil(1.5 * lambda(ldx) * ENDTIME);
    ptr = 1;
    blg = 0;
    scs = 0;
    dly = 0;
    cnt = 0;
    min_t = 0;
    status = 0;

    chnl_scs = scs;

    crp_flag = 0;

    crp_t = 0;
    crp_cnt = 0;
    crp_scs = 0;
    crp_idx = 0;
    crp_num = 0;

    pkt_list = zeros(num,3);
    pkt_list(:,1) = cumsum(exprnd(1/lambda(ldx),num,1));
    pkt_list(:,2) = pkt_list(:,1);

    while min_t <= ENDTIME
        if blg == 0
            pkt_list(ptr,1) = pkt_list(ptr,1) + BEB(pkt_list(ptr,3));
            pkt_list(ptr,3) = 1;
            blg = 1;
        end
        min_t_temp = pkt_list(ptr,1) + 1;
        new_pkt = sum(pkt_list(ptr+blg:end,1) < min_t_temp);
        if new_pkt > 0
            bof = 1 + 2 ^ 1 * rand(new_pkt,1);
            min_t_temp = min(min(pkt_list(ptr+blg:scs+blg+new_pkt,1)+bof)+1, min_t_temp);
            new_blg = sum(pkt_list(ptr+blg:scs+blg+new_pkt,1) < min_t_temp);
%             disp TEST1
%             disp(pkt_list(ptr+blg:scs+blg+new_blg,1));
%             disp TEST2
%             disp(bof(1:new_pkt,1))
            pkt_list(ptr+blg:scs+blg+new_blg,1) = pkt_list(ptr+blg:scs+blg+new_blg,1) + bof(1:new_blg);
            pkt_list(ptr+blg:scs+blg+new_blg,3) = 1;
            blg = blg + new_blg;
            pkt_list(ptr:scs+blg,:) = sortrows(pkt_list(ptr:scs+blg,:),1);
        end

        min_t = pkt_list(ptr,1) + 1;
        sect = sum(pkt_list(ptr:scs+blg,1) < min_t);
        if sect == 1
            scs = scs + 1;
            chnl_scs = chnl_scs + 1;
            dly = dly + pkt_list(ptr,1) - pkt_list(ptr,2) + 1;
            if pkt_list(ptr,1) - pkt_list(ptr,2) < 0
                disp FALSE_DLY_1
            end
            pkt_list(ptr,3) = -1;
            ptr = ptr + 1;
            blg = blg - 1;
        else
            coll_start_t = pkt_list(ptr,1);
            blg_end = sect - 1;
            for ii = 1:blg_end
                if pkt_list(ptr+ii-1,1) > pkt_list(ptr+ii,1)
                    disp FALSE_Debug0
                end
            end
            while blg_end > 0 && min_t < ENDTIME
                min_t = pkt_list(ptr+blg_end,1) + 1;
                new_blg = sum(pkt_list(ptr+blg:end,1) < min_t);
                if new_blg > 0
                    pkt_list(ptr+blg:scs+blg+new_blg,1) = pkt_list(ptr+blg:scs+blg+new_blg,1) + 1 + 2 ^ 1 * rand(new_blg,1);
                    pkt_list(ptr+blg:scs+blg+new_blg,3) = 1;
                    blg = blg + new_blg;
                    pkt_list(ptr:scs+blg,:) = sortrows(pkt_list(ptr:scs+blg,:),1);
                    for ii = 1:blg_end
                        if pkt_list(ptr+ii-1,1) > pkt_list(ptr+ii,1)
                            disp FALSE_Debug1
                        end
                    end
                end
                for ii = 1:blg_end
                    if pkt_list(ptr+ii-1,1) > pkt_list(ptr+ii,1)
                        disp FALSE_Debug2
                    end
                end
                sect = sum(pkt_list(ptr:scs+blg,1) < min_t);
                if sect == blg_end + 1
                    break;
                else
                    blg_end = sect - 1;
                end
            end

            crp_cnt = crp_cnt + 1;
            crp_start_t = min_t;
            crp_list = pkt_list(ptr:ptr+blg_end,:);
            crp_list(:,1) = min_t;
            crp_blg = sum(crp_list(:,3) >= 1);

            % F
            crp_list(1,3) = -1;
            scs = scs + 1;
            crp_scs = crp_scs + 1;
            dly = dly + crp_list(1,1) - crp_list(1,2) + 1;
            crp_list(end,1) = crp_list(1,1) + 1;

            % L
            crp_list(end,3) = -1;
            scs = scs + 1;
            crp_scs = crp_scs + 1;
            dly = dly + crp_list(end,1) - crp_list(end,2) + 1;
            crp_min_t = crp_list(end,1) + 1;

            if crp_list(end,1) - crp_list(end,2) < 0 || crp_list(1,1) - crp_list(1,2) < 0
                disp FALSE_DLY_2
            end

            if (ceil(crp_min_t) - ceil(crp_start_t)) ~= 2
                disp FALSE_SCHEDULING
            end
            crp_t = crp_t + crp_min_t - crp_start_t;

            pkt_list(ptr,:) = crp_list(1,:);
            pkt_list(ptr+1,:) = crp_list(end,:);

            % O
            if crp_blg > 2
                crp_list(2:end-1,1) = crp_list(2:end-1,1) + BEB(crp_list(2:end-1,3));
                crp_list(2:end-1,3) = crp_list(2:end-1,3) + 1;
                pkt_list(ptr+2:ptr+blg_end,:) = crp_list(2:end-1,:);
                pkt_list(ptr+2:ptr+1+blg,:) = sortrows(pkt_list(ptr+2:ptr+1+blg,:),1);
            end

            ptr = scs + 1;
            blg = blg - 2;

            prev_min_t = min_t;
            while sum(pkt_list(ptr:end,1) < crp_min_t) > 0
                if blg == 0
                    pkt_list(ptr,1) = pkt_list(ptr,1) + BEB(pkt_list(ptr,3));
                    pkt_list(ptr,3) = 1;
                    blg = 1;
                end
                min_t_temp = pkt_list(ptr,1) + 1;
                new_pkt = sum(pkt_list(ptr+blg:end,1) < min_t_temp);
                if new_pkt > 0
                    bof = 1 + 2 ^ 1 * rand(new_pkt,1);
                    min_t_temp = min(min(pkt_list(ptr+blg:scs+blg+new_pkt,1)+bof)+1, min_t_temp);
                    new_blg = sum(pkt_list(ptr+blg:scs+blg+new_pkt,1) < min_t_temp);
                    pkt_list(ptr+blg:scs+blg+new_blg,1) = pkt_list(ptr+blg:scs+blg+new_blg,1) + bof(1:new_blg);
                    pkt_list(ptr+blg:scs+blg+new_blg,3) = 1;
                    blg = blg + new_blg;
                    pkt_list(ptr:scs+blg,:) = sortrows(pkt_list(ptr:scs+blg,:),1);
                end

                min_t = pkt_list(ptr,1) + 1;
                if min_t > crp_min_t
                    min_t = prev_min_t;
                    break;
                end
                sect = sum(pkt_list(ptr:scs+blg,1) < min_t);
                if sect == 1
                    scs = scs + 1;
                    chnl_scs = chnl_scs + 1;
                    dly = dly + pkt_list(ptr,1) - pkt_list(ptr,2) + 1;
                    if pkt_list(ptr,1) - pkt_list(ptr,2) < 0
                        disp FALSE_DLY_3
                    end
                    pkt_list(ptr,3) = -1;
                    ptr = ptr + 1;
                    blg = blg - 1;
                else
                    coll_start_t = pkt_list(ptr,1);
                    blg_end = sect - 1;
                    while blg_end > 0 && min_t < ENDTIME
                        min_t = pkt_list(ptr+blg_end,1) + 1;
                        new_blg = sum(pkt_list(ptr+blg:end,1) < min_t);
                        if new_blg > 0
                            pkt_list(ptr+blg:scs+blg+new_blg,1) = pkt_list(ptr+blg:scs+blg+new_blg,1) + 1 + 2 ^ 1 * rand(new_blg,1);
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
                    if min_t > crp_min_t
                        min_t = prev_min_t;
                        break;
                    else
                        min_t = prev_min_t;
                        pkt_list(ptr:ptr+blg_end,1) = min_t + BEB(pkt_list(ptr:ptr+blg_end,3));
                        pkt_list(ptr:ptr+blg_end,3) = pkt_list(ptr:ptr+blg_end,3) + 1;
                        pkt_list(ptr:scs+blg,:) = sortrows(pkt_list(ptr:scs+blg,:),1);
                    end
                end
            end

        end
    end
    if chnl_scs + crp_scs ~= scs
        disp(ldx);
    end
    sys_thrpt(ldx) = (scs / min_t) / 2;
    chanl_thrpt(ldx) = chnl_scs / min_t;
    dly_list(ldx) = dly / scs;
    crp_thrpt(ldx) = crp_scs / min_t;
    crp_chnl_util(ldx) = crp_t / min_t;
end
toc

pt = find(sys_thrpt == max(sys_thrpt));
disp(sys_thrpt(pt));

lambda = lambda / 2;

figure
plot(lambda,sys_thrpt,lambda,chanl_thrpt,lambda,crp_thrpt,'LineWidth',1.5)
legend('System Thrpughput','Channel Throughput','CRP Throughput','Location','southeast','Interpreter','latex','FontSize',14.4)
grid on
% xlim([0 0.36])
xlabel('$\lambda$','Interpreter','latex','FontSize',17.6)
ylabel('Throughput (packet/sec)','Interpreter','latex','FontSize',17.6)
title('CRP','Interpreter','latex','FontSize',17.6)

figure
plot(lambda,dly_list,'LineWidth',1)
legend('Delay','Location','northwest','Interpreter','latex','FontSize',14.4)
grid on
ylim([0 300])
xlabel('$\lambda$','Interpreter','latex','FontSize',17.6)
ylabel('Delay (sec)','Interpreter','latex','FontSize',17.6)
title('CRP','Interpreter','latex','FontSize',17.6)

figure
plot(lambda,crp_chnl_util,'LineWidth',1)
legend('Utilization','Location','northwest','Interpreter','latex','FontSize',14.4)
grid on
xlabel('$\lambda$','Interpreter','latex','FontSize',17.6)
ylabel('Utilization','Interpreter','latex','FontSize',17.6)
title('CRP','Interpreter','latex','FontSize',17.6)

function backoff = BEB(sz_cnt)
%BEB - Description
%
% Syntax: backoff = BEB(sz_cnt)
%
% Long description
    backoff = sz_cnt;
    sz_cnt = sz_cnt + 1;
    for idx = 1:size(sz_cnt)
        backoff(idx) = 1 + 2 ^ min(sz_cnt(idx),8) * rand;
    end
end