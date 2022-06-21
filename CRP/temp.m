% New arrival in channel
while pkt_list(ptr,1) < crp_min_t
    min_t = pkt_list(ptr,1) + 1;    % packet length equals 1
    new_pkt = sum(pkt_list(ptr+blg:end,1) < min_t);
    if new_pkt > 0
        bof = exprnd(1/mu, new_pkt, 1);
        min_t = min(min(pkt_list(ptr+blg:scs+blg+new_pkt,1)+bof)+1, min_t);
        new_blg = sum(pkt_list(ptr+blg:scs+blg+new_pkt,1) < min_t);
        pkt_list(ptr+blg:scs+blg+new_blg,1) = pkt_list(ptr+blg:scs+blg+new_blg,1) + bof(1:new_blg);
        pkt_list(ptr+blg:scs+blg+new_blg,3) = 1;
        blg = blg + new_blg;
        pkt_list(ptr:scs+blg,:) = sortrows(pkt_list(ptr:scs+blg,:),1);
    end
    min_t = pkt_list(ptr,1) + 1;
    sect = sum(pkt_list(ptr:ptr+blg-1,1) < min_t);
    if sect == 1
        scs = scs + 1;
        dly = dly + pkt_list(ptr,1) - pkt_list(ptr,2) + 1;
        pkt_list(ptr,3) = -1;   % column 3 == -1 => scs
        ptr = ptr + 1;
        blg = blg - 1;
    else
        blg_end = sect - 1;

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
        pkt_list(ptr:ptr+blg_end,1) = pkt_list(ptr:ptr+blg_end,1) + exprnd(1/mu,sect,1);
        pkt_list(ptr:scs+blg,:) = sortrows(pkt_list(ptr:scs+blg,:),1);
    end
    if sum(pkt_list(:,2) < min_t) - scs ~= blg
        disp FALSE_BLG2
        return
    end
end