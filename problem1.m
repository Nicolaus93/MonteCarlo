%% initialize parameters
clear

load('stations.mat')
global Z phi psiZ psiW trans stations
deltaT = .5;
alpha = 0.6;
definePars(deltaT, alpha)

%% Problem 1: get a trajectory

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

for i = 3:500
    index = randsample(5,1,true,trans(index,:));
    X(:,i) = phi*X(:,i-1) + psiZ*Z(:,index) + ...
        psiW * mvnrnd(zeros(2,1),.5^2*ones(2))';
end

x = X(1,:);
y = X(4,:);
hold on
scatter(x,y,1,'b')
scatter(stations(1,:),stations(2,:),'r','filled')
%xlabel('X1')
%ylabel('X2')

% generate observation 
for i = 1:length(X)
    observ(:,i) = obs(x(i), y(i));  
end

%% a test on our trajectory (SISR)

N = 10000;
[tau, ~] = fastSISR(N, observ);
x1 = tau(1,:);
y1 = tau(4,:);
hold on
scatter(stations(1,:),stations(2,:),'r')
plot(x,y,'b')
plot(x1,y1,'r')

%% problem 3
load('RSSI-measurements.mat')
[tau, w] = fastSIS(Y);
x1 = tau(1,:);
y1 = tau(4,:);

figure(1)
hold on
scatter(stations(1,:),stations(2,:),'r','filled')
plot(x1,y1,'-.r')
legend('station 1', 'estimated traj','true traj')

figure(2)
l = linspace(1,60,60);
ess = effSampleSize(w(:,1:60));
scatter(l,ess,'b','filled')

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
N = 10000;
[tau, w] = fastSISR(N, Y);
x1 = tau(1,:);
y1 = tau(4,:);

% plot results
figure(1)
hold on
scatter(stations(1,:),stations(2,:),'r','*')
scatter(x1,y1,1,'b','filled')
legend('stations','estimated trajectory')

% eff sample size
figure(2)
plot(effSampleSize(w))

%% problem5 with our trajectory

c = zeros(1,10);
est = zeros(1,10);
N = 10000;

%% run this part changing alpha and c(i) values
tic
alpha = .6; % change at every iteration
phiTilde = [1, deltaT, deltaT^2/2; 0, 1, deltaT; 0, 0, alpha];
phi = [phiTilde, zeros(3); zeros(3), phiTilde];
for i = 1:10
    [~, w] = fastSISR(N, observ);          
    est(i) = sum(log(sum(w)/N))/length(observ);
end
c(6) = sum(est)/length(est); % change at every iteration
toc 

%% it does everything with one for loop, takes too much time
load('RSSI-measurements-unknown-alpha.mat')
cio = zeros(1,9);
k = 10; % num of iterations for each alpha
N = 100;

for i = 1:9
    alpha = i/10;
    phiTilde = [1, deltaT, deltaT^2/2; 0, 1, deltaT; 0, 0, alpha];
    phi = [phiTilde, zeros(3); zeros(3), phiTilde];
    est = zeros(1,k);
    for j = 1:k        
        [~, w] = fastSISR(N, Y);          
        est(j) = sum(log(sum(w)/N))/length(Y);
    end
    cio(i) = sum(est)/length(est);
end 

[~, index] = max(cio);
realAlpha = index/10
