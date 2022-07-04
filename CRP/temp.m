clear

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