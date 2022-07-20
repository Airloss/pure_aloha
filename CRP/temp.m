pkt_list(ptr(kk):ptr(kk)+blg_end(kk),1,kk) = min_t(kk) + exprnd(1/mu(kk),blg_end(kk)+1,1);
pkt_list(ptr(kk):scs(kk)+blg(kk),:,kk) = sortrows(pkt_list(ptr(kk):scs(kk)+blg(kk),:,kk),1);