%%global phi psiZ psiW
global Z phi psiZ psiW trans
deltaT = .5;
alpha = .6;
definePars(deltaT, alpha)

%%

% initial X, and W
mu = zeros(6,1);
sigma= diag([500,5,5,200,5,5]);
X(:,1) = mvnrnd(mu,sigma);
W = mvnrnd(zeros(2,1),.5^2*ones(2));

% state with n = 0 (X1)
index = randi(5);
X(:,2) = phi*X(:,1) + psiZ*Z(:,index) + ...
        psiW * mvnrnd(zeros(2,1),.5^2*ones(2))';
temp = zeros(1,5);
temp(randi(5)) = 1;

%index = randsample(5,1,true,trans(:,randi(5)))

for i = 3:100
    index = randsample(5,1,true,temp*trans);
    temp = zeros(1,5);
    temp(index) = 1;
    X(:,i) = phi*X(:,i-1) + psiZ*Z(:,index) + ...
        psiW * mvnrnd(zeros(2,1),.5^2*ones(2))';
end

x = X(1,:);
y = X(4,:);
scatter(x,y,'b')

%% observation model

load('stations.mat')
% observations parameters
for i = 1:100
    observ(:,i) = obs(x(i), y(i), stations);  
end

%%

[tau, w] = SIS(observ, stations);
x1 = tau(1,:);
y1 = tau(4,:);
hold on
scatter(x,y,'b')
scatter(x1,y1,'r')

%%
load('RSSI-measurements.mat')

