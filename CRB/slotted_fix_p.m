clear

prob = 0.1;
N = 100;

lambda = 0.01:0.01:0.4;

thr_list = zeros(length(lambda),1);
dly_list = zeros(length(lambda),1);

tic
for ldx = 1:length(lambda)
    scs = 0;
    dly = 0;
    blg = 10;
    mu = 1 / blg;
    
end
toc