kappa = 3.4;
sigma =0.59;
alpha = 0.803;
beta = 0.05;

% syms x

Xmax= 1.3;
% Xmin=double(solve(1/2*sigma^2*exp(x)+kappa*(alpha-x-1)*exp(x) - beta*(exp(x))==0,x))-1;
Xmin = 0.3;

Qmax=100;
Qmin=0;

NumX = 21;
NumQ = 21;

XVec = linspace(Xmin,Xmax,NumX);
QVec = linspace(Qmin,Qmax,NumQ);

k = 20;

plot(QVec,result{7}{k}(:,1))

figure

plot(QVec,result{7}{k}(:,2))

figure
plot(QVec,result{7}{k}(:,9))

figure
plot(QVec,result{7}{k}(:,10))

