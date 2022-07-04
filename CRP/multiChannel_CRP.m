clear

THEATA = 0.99;
ENDTIME = 1e5;
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
    num = ceil(1.5 * lambda(ldx) * ENDTIME);
    ptr = ones(CHANNEL,1);
    blg = zeros(CHANNEL,1);
    scs = blg;
    dly = blg;
    cnt = blg;
    min_t = blg;
    mu = 0.6468 / 10;

    pkt_list = zeros(num,3,CHANNEL);
    for idx = 1:CHANNEL
        pkt_list(:,1,idx) = cumsum(exprnd(1/lambda(ldx), num, 1));
        pkt_list(:,2,idx) = pkt_list(:,1,idx);
    end
    
    while min_t(1) < ENDTIME || min_t(2) < ENDTIME
        for idx = 1:2
            if blg(idx) == 0
                pkt_list(ptr(idx),1,idx) = pkt_list(ptr(idx),1,idx) + exprnd(mu,1);
                pkt_list(ptr(idx),3,idx) = 1;    % stack in backlog list
                blg(idx) = 1;
            end
            min_t_temp = pkt_list(ptr(idx),1,idx) + 1;    % packet length equals 1
            new_pkt = sum(pkt_list(ptr(idx)+blg(idx):end,1,idx) < min_t_temp);
            if new_pkt > 0
                bof = exprnd(1/mu, new_pkt, 1);
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
        end
        for idx = 1:2
            sect = sum(pkt_list(ptr(idx):scs(idx)+blg(idx),1,idx) < min_t(idx));
            if sect == 1
                scs(idx) = scs(idx) + 1;
                dly(idx) = dly(idx) + pkt_list(ptr(idx),1,idx) - pkt_list(ptr(idx),2,idx) + 1;
                pkt_list(ptr(idx),3,idx) = -1;   % column 3 == -1 => scs(idx)
                ptr(idx) = ptr(idx) + 1;
                blg(idx) = blg(idx) - 1;
            else
                coll_start_t = pkt_list(ptr(idx),1,idx);
                blg_end = sect - 1;
                cnt_coll = cnt_coll + 1;
                % collision period
                while blg_end > 0 && min_t(idx) < ENDTIME
                    min_t(idx) = pkt_list(ptr(idx)+blg_end,1,idx) + 1;
                    new_blg = sum(pkt_list(ptr(idx)+blg(idx):end,1,idx) < min_t(idx));
                    if new_blg > 0
                        pkt_list(ptr(idx)+blg(idx):scs(idx)+blg(idx)+new_blg,1,idx) = pkt_list( ...
                            ptr(idx)+blg(idx):scs(idx)+blg(idx)+new_blg,1,idx) + exprnd(1/mu,new_blg,1);
                        pkt_list(ptr(idx)+blg(idx):scs(idx)+blg(idx)+new_blg,3,idx) = 1;
                        blg(idx) = blg(idx) + new_blg;
                        pkt_list(ptr(idx):scs(idx)+blg(idx),:,idx) = sortrows(pkt_list(ptr(idx):scs(idx)+blg(idx),:,idx),1);
                    end
                    sect = sum(pkt_list(ptr(idx):scs(idx)+blg(idx),1,idx) < min_t(idx));
                    if sect == blg_end + 1
                        break;
                    else
                        blg_end = sect - 1;
                    end
                end
                
            end
        end
    end
end
toc