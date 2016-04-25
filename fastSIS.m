function [tau, w] = fastSIS(obs, resampling)
    global Z phi psiZ psiW trans 
    N = 10000;
    n = length(obs);
    tau = zeros(6,n);
    n = n-1;
    
    % initial particle
    part = mvnrnd(zeros(6,1), diag([500,5,5,200,5,5]), N)'; % 6XN particles
    
    % initial weights
    w(:,1) = prob(obs(:,1)', part(1,:), part(4,:));
    
    indices = randi(5,1,N);
    % should define tau(1)
    % tau(:,1) = 
    bigZ = zeros(2,N);
    
    % nice version
    for k = 1:n
        for j = 1:5
            temp = indices;
            indices(indices==j) = randsample(5,sum(indices==j),true,trans(j,:));
            bigZ(:,temp==j) = Z(:,indices(temp==j));
        end
        part = phi*part + psiZ*bigZ + psiW * mvnrnd(zeros(2,1),.5^2*ones(2),N)'; % mutation
        w(:,k+1) = w(:,k).*prob(obs(:,k+1)', part(1,:), part(4,:)); % weighting
        k %print the timesteps
        tau(:,k+1) = sum(bsxfun(@times,part,w(:,k+1)'),2)/sum(w(:,k+1)); % estimation
        if resampling == true
           % to do: modify w update as well
           ind = randsample(N,N,true,w(:,k+1));
           part = part(:,ind);
        end
    end