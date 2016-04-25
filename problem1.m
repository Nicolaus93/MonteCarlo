%% initialize parameters
global Z phi psiZ psiW trans
deltaT = .5;
alpha = .6;
definePars(deltaT, alpha)

%% get a trajectory

% X0 state
mu = zeros(6,1);
sigma= diag([500,5,5,200,5,5]);
X(:,1) = mvnrnd(mu,sigma);
%W = mvnrnd(zeros(2,1),.5^2*ones(2));

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
scatter(x,y,'b')

%% generate observation 

load('stations.mat')
global stations
for i = 1:length(X)
    observ(:,i) = obs(x(i), y(i), stations);  
end
hold on
scatter(stations(1,:),stations(2,:),'r')
plot(x,y,'b')

%% a test on our trajectory
tic
[tau, w] = fastSISR(observ);
toc
x1 = tau(1,:);
y1 = tau(4,:);
hold on
scatter(stations(1,:),stations(2,:),'r')
plot(x,y,'b')
plot(x1,y1,'r')

%% plot histograms
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
load('RSSI-measurements.mat')
[tau, w] = fastSISR(Y);
x1 = tau(1,:);
y1 = tau(4,:);

%% plot results
hold on
scatter(stations(1,:),stations(2,:),'r','filled')
scatter(x1,y1,1,'b','filled')
