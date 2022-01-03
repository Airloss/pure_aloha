clear

THEATA = 0.99;
ENDTIME = 1e5;

lambda = 0.01:0.01:0.4;

thrpt_list = zeros(length(lambda),1);
dly_list = zeros(length(lambda),1);

for ldx = 1:length(lambda)
    num = ceil(2.5 * lambda(ldx) * ENDTIME);
    ptr = 1;
    blg = 0;
    scs = 0;
    dly = 0;
    cnt = 0;
    min_t = 0;
    es_blg = 10;
    ac_blg = 10;
    lambda_recur = 0;
    prev_end_t = 0;
    mu = 0.83 / es_blg;

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
        min_t = pkt_list(ptr,1) + 1;

        % cycle period
        sect = pkt_list(ptr:scs+blg,1) < min_t;
        if sum(sect) == 1
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
            blg_end = sum(sect) - 1;

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
            coll_t = min_t - coll_start_t;

            % CRP period
            pkt_list(ptr:ptr+blg_end,1) = min_t;
            crp_start_t = min_t;
            crp_list = pkt_list(ptr:ptr+blg_end,:);
            if sum(crp_list(:,3) == 1) <= 1
                disp FALSE_COLL
            end
            crp_t = 0;
            crp_mu = 0.83 / 3;
            while sum(crp_list(:,3) == 1)  > 0
                if sum(crp_list(:,3) == 1) == 1
                    idx = find(crp_list(:,3) == 1);
                    crp_list(idx,3) = -1;   % column 3 == -1 => scs
                    scs = scs + 1;
                    dly = dly + crp_list(idx,1) - crp_list(idx,2) + 1;
                    blg = blg - 1;
                else
                    
                end
            end

        end

    end
end