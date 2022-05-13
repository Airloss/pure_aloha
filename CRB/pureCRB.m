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
    blg = 10;
    mu = 1 / blg;
    lambda_recur = 0.1;
    arrival_list = poissrnd(lambda(ldx) / N,N,ENDTIME);
    trans_list = zeros(N,3);    % transmit time | arrival time | synchronize state
    queue = zeros(N,ceil(1.5 * lambda(ldx) * ENDTIME));

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

        
    end
end
toc