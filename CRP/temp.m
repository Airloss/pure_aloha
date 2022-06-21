min_t_temp = min(min(pkt_list(ptr_+blg_:ptr_+blg_-1+new_pkt,1)+bof)+1, min_t_temp);
new_blg = sum(pkt_list(ptr_+blg_:ptr_+blg_-1+new_pkt,1) < min_t_temp);
pkt_list(ptr_+blg_:ptr_+blg_-1+new_blg,1) = pkt_list(ptr_+blg_:ptr_+blg_-1+new_blg,1) + bof(1:new_blg);
pkt_list(ptr_+blg_:ptr_+blg_-1+new_blg,3) = 1;
blg_ = blg_ + new_blg;
pkt_list(ptr_:ptr_+blg_-1,:) = sortrows(pkt_list(ptr_:ptr_+blg_-1,:),1);