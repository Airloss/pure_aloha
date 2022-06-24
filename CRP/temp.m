clear

channel_flag = zeros(2,1);
flag = flag_status(channel_flag);

function flag = flag_status(channel_flag)
%flag_status - Description
%
% Syntax: flag = flag_status(channel_flag)
%
% Long description
    num = sum(channel_flag == 0);
    if num == 2
        flag = 1;
        channel_flag(flag) = 1;
    elseif num == 1
        flag = find(channel_flag == 0);
        channel_flag(flag) = 1;
    else
        flag = 0;
    end
end