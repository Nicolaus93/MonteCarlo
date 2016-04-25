%% initialize parameters
load('stations.mat')
global Z phi psiZ psiW trans stations
deltaT = .5;
alpha = .6;
definePars(deltaT, alpha)

%% get a trajectory

% X0 state
mu = zeros(6,1);
sigma= diag([500,5,5,200,5,5]);
X(:,1) = mvnrnd(mu,sigma);

% state with n = 0 (X1)
index = randi(5);
X(:,2) = phi*X(:,1) + psiZ*Z(:,index) + ...
        psiW * mvnrnd(zeros(2,1),.5^2*ones(2))';
temp = zeros(1,5);
temp(randi(5)) = 1;

%index = randsample(5,1,true,trans(:,randi(5)))

for i = 3:500
    index = randsample(5,1,true,trans(index,:));
    %index = randsample(5,1,true,temp*trans);
    %temp = zeros(1,5);
    %temp(index) = 1;
    X(:,i) = phi*X(:,i-1) + psiZ*Z(:,index) + ...
        psiW * mvnrnd(zeros(2,1),.5^2*ones(2))';
end

x = X(1,:);
y = X(4,:);
hold on
scatter(x,y,1,'b')
scatter(stations(1,:),stations(2,:),'r','filled')


%% generate observation 

for i = 1:length(X)
    observ(:,i) = obs(x(i), y(i));  
end

%% a test on our trajectory (SISR)

[tau, ~] = fastSISR(observ);
x1 = tau(1,:);
y1 = tau(4,:);
hold on
scatter(stations(1,:),stations(2,:),'r')
plot(x,y,'b')
plot(x1,y1,'r')

%% problem 3
load('RSSI-measurements.mat')
[tau, w] = fastSIS(Y, false);
x1 = tau(1,:);
y1 = tau(4,:);
hold on
scatter(stations(1,:),stations(2,:),'r')
plot(x1,y1,'r')

%% plot histograms for weights
figure
subplot(4,1,1)       
histogram(log10(w(:,1)),-350:10:0)
title('n = 1')

subplot(4,1,2)       
histogram(log10(w(:,10)),-350:10:0)
title('n = 10')

subplot(4,1,3)       
histogram(log10(w(:,20)),-350:10:0)
title('n = 20')

subplot(4,1,4)       
histogram(log10(w(:,40)),-350:10:0)
title('n = 40')

%% problem4 trajectory estimation
[tau, w] = fastSISR(Y);
x1 = tau(1,:);
y1 = tau(4,:);

%% plot results
hold on
scatter(stations(1,:),stations(2,:),'r','filled')
scatter(x1,y1,1,'b','filled')
