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
    % index = randsample(5,1,true,trans(index,:);
    index = randsample(5,1,true,temp*trans);
    temp = zeros(1,5);
    temp(index) = 1;
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

%% a test on our trajectory
%tic
%[tau1, w1] = SIS(observ, stations);
%toc
tic
[tau, w] = fastSIS(observ);
toc
x1 = tau(1,:);
y1 = tau(4,:);
hold on
scatter(x,y,'b')
plot(x1,y1,'r')

%% to do
load('RSSI-measurements.mat')

