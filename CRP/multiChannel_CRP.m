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

    pkt_list = zeros(num,3,CHANNEL);
    for idx = 1:CHANNEL
        pkt_list(:,1,idx) = cumsum(exprnd(1/lambda(ldx), num, 1));
        pkt_list(:,2,idx) = pkt_list(:,1,idx);
    end
    
    while min_t(1) < ENDTIME || min_t(2) < ENDTIME
        for idx = 1:2
            if blg(idx) == 0
                pkt_list(ptr,1,idx) = pkt_list(ptr,1,idx) + exprnd(mu,1,idx);
                pkt_list(ptr,3,idx) = 1;    % stack in backlog list
                blg(idx) = 1;
            end
            min_t_temp = pkt_list(ptr,1,idx) + 1;    % packet length equals 1
            new_pkt = sum(pkt_list(ptr+blg:end,1,idx) < min_t_temp);
            if new_pkt > 0
                bof = exprnd(1/mu, new_pkt, 1);
                min_t_temp = min(min( ...
                    pkt_list(ptr+blg:scs+blg+new_pkt,1,idx)+bof)+1, min_t_temp);
                new_blg = sum(pkt_list(ptr+blg:scs+blg+new_pkt,1,idx) < min_t_temp);
                pkt_list(ptr+blg:scs+blg+new_blg,1,idx) = pkt_list( ...
                    ptr+blg:scs+blg+new_blg,1,idx) + bof(1:new_blg);
                pkt_list(ptr+blg:scs+blg+new_blg,3,idx) = 1;
                blg = blg + new_blg;
                pkt_list(ptr:scs+blg,:,idx) = sortrows(pkt_list(ptr:scs+blg,:,idx),1);
            end
            % idle_t = idle_t + pkt_list(ptr,1) - min_t;
            % if idle_t < 0
            %     disp FALSE_IDLE_MINUS
            %     % return
            % end
            min_t(idx) = pkt_list(ptr,1,idx) + 1;
        end
    end
end
toc