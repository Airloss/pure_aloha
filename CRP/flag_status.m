function flag = flag_status(channel_flag)
    num = sum(channel_flag == 0);
    if num == 2
        flag = 1;
    elseif num == 1
        flag = find(channel_flag == 0);
    else
        flag = 0;
    end
end