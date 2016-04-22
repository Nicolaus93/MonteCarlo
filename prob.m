function p = prob(obs, part1, part4)
    global stations
    % stations is 2x7
    % part1=part2 is 1XN, should be a column vector
    % particles (x1 x2) should be 2xN
    % obs should be a row vector 1x7
    one = bsxfun(@minus,part1(:),stations(1,:)); % Nx7
    two = bsxfun(@minus,part4(:),stations(2,:)); % Nx7
    mean = 90-30*log10(sqrt(one.^2 + two.^2));   % Nx7
    y = repmat(obs,length(part1),1); % should be Nx7
    p = prod(normpdf(y,mean,1.5^2),2);
    %p = sum(p);
    
    %p = @(y,x) mvnpdf(y,v - ...
    %    10*eta*log10(sqrt(sum(bsxfun(@minus,x',stations).^2,1))'), ...
    %                    eye(7)*1.5^2);
