clear

THEATA = 0.99;
ENDTIME = 1e5;
N = 50;

lambda = 0.02:0.02:0.4;

thr_list = zeros(length(lambda),1);
dly_list = zeros(length(lambda),1);

tic
for ldx = 1:length(lambda)
    scs = 0;
    dly = 0;
    mu = 10;
    prob = 1 / mu;
    lambda_recur = 0.5;
    bof_list = zeros(N,1);
    arrival_list = poissrnd(lambda(ldx) / N,N,ENDTIME);
    trans_list = zeros(N,4);    % transmit time | arrival time | packets number in queue | synchronize state
    queue = zeros(N,ceil(1.5 * lambda(ldx) / N * ENDTIME));

    for tdx = 1:ENDTIME
        % Queue update before transmission
        if sum(arrival_list(:,tdx)) ~= 0
            arrival = find(arrival_list(:,tdx));
            for idx = 1:size(arrival,1)
                ndx = arrival(idx);
                append_ = sum(queue(ndx,:) > 0);
                queue(ndx,append_+1:append_+arrival_list(ndx,tdx)) = tdx;
                if append_ == 0
                    trans_list(ndx,1:2) = tdx;
                end
                trans_list(ndx,3) = trans_list(ndx,3) + arrival_list(ndx,tdx);
            end
        end

        blg = sum(trans_list(:,1) > 0 && trans_list(:,1) == tdx) > 0;
        if blg > 0
            prob_ = rand(blg,1) <= prob;
            trans_idx = trans_list(:,1) > 0;
            trans_ = trans_list(trans_idx,:) .* prob_;
            if sum(trans_(:,1) > 0) == 0
                mu = max(mu - 1, 0);
                lambda_recur = 0.995 * lambda_recur;
            elseif sum(trans_(:,1) > 0) == 1
                ndx_ = find(trans_list(:,2:4) == trans_(prob,2:4));
                scs = scs + 1;
                dly = dly + trans_(prob_,1) - trans_(prob_,2) + 1;
                trans_list(ndx_,3) = trans_list(ndx_,3) - 1;
                mu = max(mu - 1, 0);
                lambda_recur = 0.995 * lambda_recur + 0.005;
                if trans_list(ndx_,3) > 0
                    trans_list(ndx_,1:2) = queue(ndx_,2);
                    queue(ndx_,2:end) = queue(ndx_,1:end-1);
                    trans_list(ndx_,4) = trans_list(ndx_,4) + 1; % synchronized
                    % allocate slot to transmit
                    bof = ceil((2 ^ trans_list(ndx_,3)) * mu * rand);
                    while sum(bof_list == trans_list(ndx_,1) + bof) > 0
                        bof = ceil((2 ^ trans_list(ndx_,3)) * mu * rand);
                    end
                    trans_list(ndx_,1) = trans_list(ndx_,1) + bof;
                    bof_list(ndx_,1) = trans_list(ndx_,1) + bof;
                end
            else
                find(trans_(:,4) > 0 && trans_(:,1) == tdx)
            end
        end
    end
end
toc