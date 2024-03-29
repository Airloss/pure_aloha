clear

THEATA = 0.99;
ENDTIME = 1e5;
SLOT = 1;

lambda = 0.01:0.01:0.30;

thrpt_list = zeros(length(lambda),1);
dly_list = zeros(length(lambda),1);
scs_crp_list = zeros(length(lambda),1); % Record the mean of success packets in CRP

tic
parfor (ldx = 1:length(lambda),6)
    num = ceil(3.5 * lambda(ldx) * ENDTIME);
    ptr = 1;
    blg = 0;
    scs = 0;
    dly = 0;
    cnt = 0;
    min_t = 0;

    pkt_list = zeros(num,3);
    pkt_list(:,1) = cumsum(exprnd(1/lambda(ldx), num, 1));
    pkt_list(:,2) = pkt_list(:,1);

    while min_t < ENDTIME
        cnt = cnt + 1;

        if blg == 0
            pkt_list(ptr,3) = 1;    % column 3 == 1 stack in backlog list
            blg = 1;
        end

        min_t = pkt_list(ptr,1) + 1;    % packet length equals 1
        new_blg = sum(pkt_list(ptr+blg:end,1) < min_t);
        if new_blg > 0
            pkt_list(ptr+blg:scs+blg+new_blg,3) = 1;
            blg = blg + new_blg;
            pkt_list(ptr:scs+blg,:) = sortrows(pkt_list(ptr:scs+blg,:),1);
        end

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

            % collision period
            while blg_end > 0 && min_t < ENDTIME
                min_t = pkt_list(ptr+blg_end,1) + 1;
                new_blg = sum(pkt_list(ptr+blg:end,1) < min_t);
                if new_blg > 0
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
            pkt_list(ptr:ptr+blg_end,3) = pkt_list(ptr:ptr+blg_end,3) + 1;  % column 3 += 1 => retransmit

            % ? collision process
            pkt_list(ptr,1) = min_t;
            pkt_list(ptr+1:ptr+blg_end-1,1) = pkt_list(ptr+blg_end,1) + 3 * SLOT + 1 + 2 .^ min((pkt_list(ptr+1:ptr+blg_end-1,3)),10) .* rand(blg_end-1,1);
            pkt_list(ptr+blg_end,1) = min_t + 1 * SLOT;
            pkt_list(ptr:scs+blg,:) = sortrows(pkt_list(ptr:scs+blg,:),1);

            if pkt_list(ptr,1) ~= min_t && pkt_list(ptr+1,1) ~= min_t + 1
                disp FALSE COLL_END
            end

            pkt_first = pkt_list(ptr,1);
            pkt_last = pkt_list(ptr+1,1);
            
            % first period
            min_t = pkt_list(ptr,1) + 1;
            new_blg = sum(pkt_list(ptr+blg:end,1) < min_t);
            if new_blg > 0
                pkt_list(ptr+blg:scs+blg+new_blg,3) = 1;
                blg = blg + new_blg;
                pkt_list(ptr:scs+blg,:) = sortrows(pkt_list(ptr:scs+blg,:),1);
            end
            sect = sum(pkt_list(ptr:scs+blg,1) < min_t);
            if sect == 1
                scs = scs + 1;
                dly = dly + pkt_list(ptr,1) - pkt_list(ptr,2) + 1;
                pkt_list(ptr,3) = -1;   % column 3 == -1 => scs
                ptr = ptr + 1;
                blg = blg - 1;

                % last period
                if pkt_list(ptr,1) == pkt_last
                    min_t = pkt_list(ptr,1) + 1;
                    new_blg = sum(pkt_list(ptr+blg:end,1) < min_t);
                    if new_blg > 0
                        pkt_list(ptr+blg:scs+blg+new_blg,3) = 1;
                        blg = blg + new_blg;
                        pkt_list(ptr:scs+blg,:) = sortrows(pkt_list(ptr:scs+blg,:),1);
                    end
                    sect = sum(pkt_list(ptr:scs+blg,1) < min_t);
                    if sect == 1
                        scs = scs + 1;
                        dly = dly + pkt_list(ptr,1) - pkt_list(ptr,2) + 1;
                        pkt_list(ptr,3) = -1;   % column 3 == -1 => scs
                        ptr = ptr + 1;
                        blg = blg - 1;
                    end
                end
            end
        end
        % pkt_drop = pkt_list(ptr:scs+blg,3) > 32;
        % if sum(pkt_drop) > 0
        %     pkt_list(pkt_drop,1) = min_t - 1 - rand(sum(pkt_drop),1);
        %     pkt_list(pkt_drop,3) = -2;
        %     pkt_list(ptr:scs+blg,:) = sortrows(pkt_list(ptr:scs+blg,:),1);
        %     blg = blg - sum(pkt_drop);
        %     scs = scs + sum(pkt_drop);
        %     ptr = scs + 1;
        % end
    end
    thrpt_list(ldx) = scs / min_t;
    dly_list(ldx) = dly / scs;
end
toc

figure
plot(lambda,thrpt_list,'LineWidth',1)
legend('Estimated','Location','northwest','Interpreter','latex','FontSize',14.4)
grid on
xlim([0 0.36])
xlabel('$\lambda$','Interpreter','latex','FontSize',17.6)
ylabel('Throughput (packet/sec)','Interpreter','latex','FontSize',17.6)
title('Pure ALOHA BEB','Interpreter','latex','FontSize',17.6)

figure
plot(lambda,dly_list,'LineWidth',1)
legend('Estimated','Location','northwest','Interpreter','latex','FontSize',14.4)
grid on
xlim([0 0.3])
ylim([0 300])
xlabel('$\lambda$','Interpreter','latex','FontSize',17.6)
ylabel('Delay (sec)','Interpreter','latex','FontSize',17.6)
title('Pure ALOHA BEB','Interpreter','latex','FontSize',17.6)