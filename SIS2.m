function [tau, w] = SIS2(obs)
    global Z phi psiZ psiW trans 
    N = 10000;
    n = length(obs);
    tau = zeros(6,n);
    n = n-1;
    
    % initial particle
    part = mvnrnd(zeros(6,1), diag([500,5,5,200,5,5]), N)'; % 6XN particles
    
    % initial weights
    w = prob(obs(:,1)', part(1,:), part(4,:));
    
    % transition matrices
    temp = zeros(N,5);
    indices = randi(5,1,N);
    for i = 1:N
        temp(i, indices(i))=1;
    end
        
    % we should get rid of one for loop here 
    for k = 1:n
        for j = 1:N
            index = randsample(5,1,true,temp(j,:)*trans);
            temp(j,:) = zeros(1,5);
            temp(j,index) = 1;
            part(:,j) = phi*part(:,j) + psiZ*Z(:,index) + ...
                psiW * mvnrnd(zeros(2,1),.5^2*ones(2))';            
        end 
        w = w.*prob(obs(:,k+1)', part(1,:), part(4,:));
        k %print the timesteps
        tau(:,k+1) = sum(bsxfun(@times,part,w'),2)/sum(w);
    end

