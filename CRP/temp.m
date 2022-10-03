clear

ENDTIME = 3e5;

lambda = 0.01:0.01:0.24;
betaT = 0.5;

thrpt_list = zeros(length(lambda),1);
dly_list = zeros(length(lambda),1);
cnt_s_list = dly_list;
cnt_c_list = dly_list;
cnt_cc_list = dly_list;

tic
parfor (ldx = 1:length(lambda),6)
    num = ceil(1.5 * lambda(ldx) * ENDTIME);
    ptr = 1;
    blg = 0;
    scs = blg;
    dly = blg;
    cnt = blg;
    min_t = blg;
    idle_t = blg;
    status = blg;
    blg_end = blg;
    mu = 0.6468 / 2;

    cnt_cc = zeros(2,1);

    pkt_list = zeros(num,3);
    pkt_list(:,1) = cumsum(exprnd(1/lambda(ldx), num, 1));
    pkt_list(:,2) = pkt_list(:,1);

    while min_t <= ENDTIME
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
            min_t_temp = min(min( ...
                pkt_list(ptr+blg:scs+blg+new_pkt,1)+bof)+1, min_t_temp);
            new_blg = sum(pkt_list(ptr+blg:scs+blg+new_pkt,1) < min_t_temp);
            pkt_list(ptr+blg:scs+blg+new_blg,1) = pkt_list( ...
                ptr+blg:scs+blg+new_blg,1) + bof(1:new_blg);
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

        sect = sum(pkt_list(ptr:scs+blg,1) < min_t);
        if sect == 1
            scs = scs + 1;
            dly = dly + pkt_list(ptr,1) - pkt_list(ptr,2) + 1;
            pkt_list(ptr,3) = -1;   % column 3 == -1 => scs
            ptr = ptr + 1;
            blg = blg - 1;
            mu = betaT / max(blg,1);
        else
            if cnt - cnt_cc(1,1) == 1
                cnt_cc(2,1) = cnt_cc(2,1) + 1;
            end
            cnt_cc(1,1) = cnt;
            coll_start_t = pkt_list(ptr,1);
            blg_end = sect - 1;
            % collision period
            while blg_end > 0 && min_t < ENDTIME
                min_t = pkt_list(ptr+blg_end,1) + 1;
                new_blg = sum(pkt_list(ptr+blg:end,1) < min_t);
                if new_blg > 0
                    pkt_list(ptr+blg:scs+blg+new_blg,1) = pkt_list( ...
                        ptr+blg:scs+blg+new_blg,1) + exprnd(1/mu,new_blg,1);
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
            mu = betaT / max(blg,1);
            pkt_list(ptr:ptr+blg_end,1) = min_t + exprnd(1/mu,sect,1);
            pkt_list(ptr:scs+blg,:) = sortrows(pkt_list(ptr:scs+blg,:),1);
        end
        if sum(pkt_list(:,2) < min_t) - scs ~= blg
            disp FALSE_BLG2
            % return
        end
        ac_blg = sum(pkt_list(:,2) < min_t) - scs;
        if ac_blg ~= blg
            disp FALSE_BLG
            % return
        end
        prev_end_t = min_t;
%         mu = betaT / max(ac_blg,1);
    end
    thrpt_list(ldx) = scs / min_t;
    dly_list(ldx) = dly / scs;
    cnt_s_list(ldx) = scs /cnt;
    cnt_c_list(ldx) = (1 - cnt_s_list(ldx))^2;
    cnt_cc_list(ldx) = cnt_cc(2,1) / cnt;
end
toc

pt = find(thrpt_list == max(thrpt_list));
disp(thrpt_list(pt));
disp(lambda(pt));

figure
plot(lambda,thrpt_list,'LineWidth',1)
grid on
xlabel('$\lambda$','Interpreter','latex','FontSize',17.6)
ylabel('Throughput (packet/sec)','Interpreter','latex','FontSize',17.6)
title('Pure ALOHA','Interpreter','latex','FontSize',17.6)

figure
plot(lambda,dly_list,'LineWidth',1)
grid on
ylim([0 300])
xlabel('$\lambda$','Interpreter','latex','FontSize',17.6)
ylabel('Throughput (packet/sec)','Interpreter','latex','FontSize',17.6)
title('Pure ALOHA','Interpreter','latex','FontSize',17.6)

figure
plot(lambda,cnt_c_list,lambda,cnt_cc_list,'LineWidth',1)
grid on
xlabel('$\lambda$','Interpreter','latex','FontSize',17.6)
ylabel('Throughput (packet/sec)','Interpreter','latex','FontSize',17.6)
title('Pure ALOHA','Interpreter','latex','FontSize',17.6)