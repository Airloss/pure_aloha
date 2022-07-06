% CRP period in addtional channel
crp_cnt = crp_cnt + 1;
pkt_list(ptr(idx):ptr(idx)+blg_end(idx),1,idx) = min_t(idx);
crp_start_t = min_t(idx);
crp_list = pkt_list(ptr(idx):ptr(idx)+blg_end(idx),:,idx);
if sum(crp_list(:,3) == 1) <= 1
    disp FALSE_COLL
end
crp_min_t = 0;
crp_blg = sum(crp_list(:,3)==1);
crp_list(:,1) = crp_list(:,1) + exprnd(1/crp_mu,crp_blg,1);
crp_list = sortrows(crp_list,1);
crp_coll = crp_coll + crp_blg;
while sum(crp_list(:,3) == 1)  > 0 && crp_min_t < ENDTIME
    if sum(crp_list(:,3) == 1) == 1
        idx = find(crp_list(:,3) == 1);
        crp_list(idx,3) = -1;   % column 3 == -1 => scs
        scs = scs + 1;
        crp_scs = crp_scs + 1;
        dly = dly + crp_list(idx,1) - crp_list(idx,2) + 1;
        blg = blg - 1;
        crp_min_t = crp_list(idx,1) + 1;
    else
        crp_trans_idx = crp_list(:,3) == 1;
        crp_trans_list = crp_list(crp_trans_idx,:);
        crp_min_t = crp_trans_list(1,1,idx) + 1;
        crp_sect = sum(crp_trans_list(:,1,idx) < crp_min_t);
        if crp_sect == 1
            scs = scs + 1;
            crp_scs = crp_scs + 1;
            dly = dly + crp_trans_list(1,1,idx) - crp_trans_list(1,2) + 1;
            crp_trans_list(1,3) = -1;
            blg = blg - 1;
            crp_trans_list(2:end,1,idx) = crp_trans_list(1,1,idx) + 1;
        else
            crp_sect_num = crp_sect - 1;
            while crp_sect_num > 0 && crp_min_t < ENDTIME
                crp_min_t = crp_trans_list(crp_sect,1,idx) + 1;
                crp_sect = sum(crp_trans_list(:,1,idx) < crp_min_t);
                crp_trans_list(1:crp_sect,3) = 2;
                if crp_sect == crp_sect_num + 1
                    break;
                else
                    crp_sect_num = crp_sect - 1;
                end
            end
            crp_trans_list(1:crp_sect,1,idx) = crp_min_t ...
                + exprnd(1/crp_mu,crp_sect,1);
            crp_trans_list(:,3) = crp_trans_list(:,3) - 1;  % particapate: 1, not: 0
            crp_trans_list = sortrows(crp_trans_list,1,idx);
        end
        crp_list(crp_trans_idx,:) = crp_trans_list;
    end
end
crp_t = crp_t + crp_min_t - crp_start_t;

if sum(crp_list(:,3) < -1) > 0
    disp FALSE_CRP_SCS
end