function y = obs(x,y)
    global stations
    v = 90*ones(length(stations),1);
    eta = 3;
    diff = bsxfun(@minus,[x,y]',stations);
    eucDist = sqrt(sum(diff.^2,1))';
    y = v - 10*eta*log10(eucDist) + ... 
        + mvnrnd(zeros(length(stations),1),...
        1.5*ones(length(stations)))';
    