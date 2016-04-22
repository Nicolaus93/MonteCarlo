function [tau, w] = SIS(obs,stations)
    global Z phi psiZ psiW trans 
    N = 100;
    n = length(obs)-1;
    tau = zeros(6,n);
    
    % prob dist of w
    v = 90*ones(length(stations),1);
    eta = 3;
    p = @(y,x) mvnpdf(y,v - ...
        10*eta*log10(sqrt(sum(bsxfun(@minus,x',stations).^2,1))'), ...
                        eye(7)*1.5^2);
                    
    % initial particle
    part = mvnrnd(zeros(6,1), diag([500,5,5,200,5,5]), N)'; % 6XN particles
    
    % weights
    w = zeros(N,1);
    for j = 1:N
        w(j) = p(obs(:,1), [part(1,j), part(4,j)]);
    end
    tau(:,1) = sum(bsxfun(@times,part,w'),2)/sum(w);
    
    % transition matrices
    temp = zeros(N,5);
    indices = randi(5,1,N);
    for i = 1:N
        temp(i, indices(i))=1;
    end
    
    for k = 1:n
        for j = 1:N
            index = randsample(5,1,true,temp(j,:)*trans);
            temp(j,:) = zeros(1,5);
            temp(j,index) = 1;
            part(:,j) = phi*part(:,j) + psiZ*Z(:,index) + ...
                psiW * mvnrnd(zeros(2,1),.5^2*ones(2))';
            w(j) = w(j).*p(obs(:,k+1), [part(1,j), part(4,j)]);
        end  
        tau(:,k+1) = sum(bsxfun(@times,part,w'),2)/sum(w);
    end
    