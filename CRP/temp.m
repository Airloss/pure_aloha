for ii = 1:2
    while sum(pkt_list(ptr(ii):end,1,ii) < crp_min_t) > 0
        if blg(ii) == 0
            pkt_list(ptr(ii),1,ii) = pkt_list(ptr(ii),1,ii) + exprnd(mu(ii),1);
            pkt_list(ptr(ii),3,ii) = 1;    % stack in backlog list
            blg(ii) = 1;
        end
        min_t_temp = pkt_list(ptr(ii),1,ii) + 1;    % packet length equals 1
        new_pkt = sum(pkt_list(ptr(ii)+blg(ii):end,1,ii) < min_t_temp);
        if new_pkt > 0
            bof = exprnd(1/mu(ii), new_pkt, 1);
            min_t_temp = min(min( ...
                pkt_list(ptr(ii)+blg(ii):scs(ii)+blg(ii)+new_pkt,1,ii)+bof)+1, min_t_temp);
            new_blg = sum(pkt_list(ptr(ii)+blg(ii):scs(ii)+blg(ii)+new_pkt,1,ii) < min_t_temp);
            pkt_list(ptr(ii)+blg(ii):scs(ii)+blg(ii)+new_blg,1,ii) = pkt_list( ...
                ptr(ii)+blg(ii):scs(ii)+blg(ii)+new_blg,1,ii) + bof(1:new_blg);
            pkt_list(ptr(ii)+blg(ii):scs(ii)+blg(ii)+new_blg,3,ii) = 1;
            blg(ii) = blg(ii) + new_blg;
            pkt_list(ptr(ii):scs(ii)+blg(ii),:,ii) = sortrows(pkt_list(ptr(ii):scs(ii)+blg(ii),:,ii),1);
        end
        % idle_t = idle_t + pkt_list(ptr(ii),1) - min_t;
        % if idle_t < 0
        %     disp FALSE_IDLE_MINUS
        %     % return
        % end
        min_t(ii) = pkt_list(ptr(ii),1,ii) + 1;
        sect = sum(pkt_list(ptr(ii):scs(ii)+blg(ii),1,ii) < min_t(ii));
        if sect == 1
            scs(ii) = scs(ii) + 1;
            dly(ii) = dly(ii) + pkt_list(ptr(ii),1,ii) - pkt_list(ptr(ii),2,ii) + 1;
            pkt_list(ptr(ii),3,ii) = -1;   % column 3 == -1 => scs(ii)
            ptr(ii) = ptr(ii) + 1;
            blg(ii) = blg(ii) - 1;
        else
            coll_start_t = pkt_list(ptr(ii),1,ii);
            blg_end(ii) = sect - 1;
            % collision period
            while blg_end(ii) > 0 && min_t(ii) < ENDTIME
                min_t(ii) = pkt_list(ptr(ii)+blg_end(ii),1,ii) + 1;
                new_blg = sum(pkt_list(ptr(ii)+blg(ii):end,1,ii) < min_t(ii));
                if new_blg > 0
                    pkt_list(ptr(ii)+blg(ii):scs(ii)+blg(ii)+new_blg,1,ii) = pkt_list( ...
                        ptr(ii)+blg(ii):scs(ii)+blg(ii)+new_blg,1,ii) + exprnd(1/mu(ii),new_blg,1);
                    pkt_list(ptr(ii)+blg(ii):scs(ii)+blg(ii)+new_blg,3,ii) = 1;
                    blg(ii) = blg(ii) + new_blg;
                    pkt_list(ptr(ii):scs(ii)+blg(ii),:,ii) = sortrows(pkt_list(ptr(ii):scs(ii)+blg(ii),:,ii),1);
                end
                sect = sum(pkt_list(ptr(ii):scs(ii)+blg(ii),1,ii) < min_t(ii));
                if sect == blg_end(ii) + 1
                    break;
                else
                    blg_end(ii) = sect - 1;
                end
            end
            pkt_list(ptr(ii):ptr(ii)+blg_end(ii),1,ii) = min_t(ii) + exprnd(1/mu(ii),blg_end(ii)+1,1);
            pkt_list(ptr(ii):scs(ii)+blg(ii),:,ii) = sortrows(pkt_list(ptr(ii):scs(ii)+blg(ii),:,ii),1);
        end
    end
end