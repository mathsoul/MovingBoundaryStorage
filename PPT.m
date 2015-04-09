%% Presentation 
tic

kappa = 3.4;
sigma =0.59;
alpha = 0.803;

Qmax=100;
Qmin=0;

% mark1
% Xmax= 2;
% Xmin = -0.6;
% beta = 0.05;
% 
% BuyingType = 'linear';
% BuyingCostPara = [0.02,0.005,-0.1]; 
% 
% SellingType = 'linear';
% SellingCostPara = [0.02,0.005,-0.1]; 

% BuyingType = 'linear';
% BuyingCostPara = [0.1,0,-0.1];
% 

%mark2
Xmax= 2.4;
Xmin = -4;

beta = 0.5;

BuyingType = 'reciprocal';
BuyingCostPara = [0.2,0,-0.2];

SellingType = 'reciprocal';
SellingCostPara = [0.2,0,-0.2];
% 
%mark3
% Xmax= 2.2;
% Xmin = -4;
% 
% beta = 0.5;
% 
% BuyingType = 'reciprocal';
% BuyingCostPara = [0.2,0,-0.2];
% 
% SellingType = 'linear';
% SellingCostPara = [0.02,0,-0.2];

% mark4
% Xmax= 2.2;
% Xmin = -4;
% beta = 0.5;
% 
% BuyingType = 'linear';
% BuyingCostPara = [0.01,0.00,-0.1]; 
% 
% SellingType = 'linear';
% SellingCostPara = [0.02,0.05,-0.1]; 


NumX = 41;
NumQ = 41;

%% Move the boundary along X at the same time
XVec = linspace(Xmin,Xmax,NumX);
QVec = linspace(Qmin,Qmax,NumQ);

Counter = 1;

clear NewPolicy BuyingProfit SellingProfit HoldingProfit Valuefunction m

M = zeros(1,NumX);
U = zeros(1,NumX);

a0 = beta/(2*kappa);
b0 = 1/2;

for i = 1:NumX
    M(i) = hypergeom(a0,b0,kappa/sigma^2*(XVec(i)-alpha)^2);
    U(i) = mchgu(a0,b0,kappa/sigma^2*(XVec(i)-alpha)^2);
end

InitialPolicy = [zeros(1,NumX);zeros(NumQ-1,NumX-2),2*ones(NumQ-1,2)];

%% Move Together
% [NewPolicy{1},Valuefunction{1},HoldingProfit{1},BuyingProfit{1},SellingProfit{1}] = ...
%     MovingBoundaryXWhole(M,U,alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,NumX,NumQ,InitialPolicy,BuyingType,BuyingCostPara,SellingType,SellingCostPara,'Together');
% 
% 
% 
% while norm(NewPolicy{Counter}-InitialPolicy)>0
%     InitialPolicy = NewPolicy{Counter};
%     Counter = Counter +1;
%    [NewPolicy{Counter},Valuefunction{Counter},HoldingProfit{Counter},BuyingProfit{Counter},SellingProfit{Counter}] = ...
%     MovingBoundaryXWhole(M,U,alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,NumX,NumQ,InitialPolicy,BuyingType,BuyingCostPara,SellingType,SellingCostPara,'Together');    
% end


 %% Move Alternatively
% dbstop in MovingBoundaryXWhole at 150 if 'MoveIndicator == 6'
dbstop if error
MoveIndicator = 1;

[NewPolicy{1},Valuefunction{1},HoldingProfit{1},BuyingProfit{1},SellingProfit{1}] = ...
    MovingBoundaryXWhole(M,U,alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,NumX,NumQ,InitialPolicy,BuyingType,BuyingCostPara,SellingType,SellingCostPara,'Alternatively',MoveIndicator);
    
MoveIndicator = MoveIndicator + 1;
Counter = Counter + 1;
[NewPolicy{2},Valuefunction{2},HoldingProfit{2},BuyingProfit{2},SellingProfit{2}] = ...
    MovingBoundaryXWhole(M,U,alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,NumX,NumQ,NewPolicy{1},BuyingType,BuyingCostPara,SellingType,SellingCostPara,'Alternatively',MoveIndicator);

oldPolicy = InitialPolicy;

while norm(NewPolicy{Counter}-oldPolicy)>0
    MoveIndicator = MoveIndicator + 1;
    oldPolicy = NewPolicy{Counter-1};
    Counter = Counter +1;
    [NewPolicy{Counter},Valuefunction{Counter},HoldingProfit{Counter},BuyingProfit{Counter},SellingProfit{Counter}] = ...
    MovingBoundaryXWhole(M,U,alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,NumX,NumQ,NewPolicy{Counter - 1},BuyingType,BuyingCostPara,SellingType,SellingCostPara,'Alternatively',MoveIndicator);    %'Alternatively'
end


IMovingBoundaryWholeX = NewPolicy{Counter};

toc



%% Plot the trajectory of moving boundary 
figure

plotBoundary(InitialPolicy,Qmax,Qmin,Xmax,Xmin,alpha,'fill')

m(1) = getframe;

for k  = 1:Counter
    Policy = NewPolicy{k};

figure

plotBoundary(Policy,Qmax,Qmin,Xmax,Xmin,alpha,'fill')

m(k+1) = getframe;
end

% movie2avi(m(1:2),'C:\Users\msbcr563\Google Drive\Longcode\presentation\MBM1','compression','none','fps',0.2)
% 
% movie2avi(m(1:2),'C:\Users\msbcr563\Google Drive\Longcode\presentation\MBM1','compression','none','fps',0.5)

% movie2avi(m,'C:\Users\msbcr563\Google Drive\Longcode\presentation\MBM','compression','none','fps',2)

% close all

%% Plot the selling profit

[~,~,~,SP] = SolvePDEWithTables(M,U,alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,NumX,NumQ,...
    InitialPolicy,BuyingType,BuyingCostPara,SellingType,SellingCostPara);

plotBoundary(InitialPolicy,Qmax,Qmin,Xmax,Xmin,alpha,'fill')

hold on
p = plot([QVec((NumQ+1)/2),QVec((NumQ+1)/2)],[Xmin,Xmax],'k');

set(p,'LineWidth',2);
hold off

figure 
plotBoundary(NewPolicy{1},Qmax,Qmin,Xmax,Xmin,alpha,'fill')

hold on
p = plot([QVec((NumQ+1)/2),QVec((NumQ+1)/2)],[Xmin,Xmax],'k');

set(p,'LineWidth',2);
hold off


figure 
SP = SP((NumQ+1)/2,:);
p = plot(XVec,SP);
titlestr = sprintf('Selling Profit with Initial Policy at q = %g',QVec((NumQ+1)/2));
title(titlestr)
set(p,'LineWidth',3);
hold on
scatter(XVec(40),0,'MarkerEdgeColor','k','MarkerFaceColor','k');
scatter(XVec(39),SP(39),'MarkerEdgeColor','g','MarkerFaceColor','g')
hold off
ylabel('Selling Profit','fontsize',13)
xlabel('Log Price x','fontsize',13)

%% Plot the buying profit
NB = 32;
BP = BuyingProfit{14}(NB,:);


% current policy
plotBoundary(NewPolicy{13},Qmax,Qmin,Xmax,Xmin,alpha,'fill')

hold on
p = plot([QVec(NB),QVec(NB)],[Xmin,Xmax],'k');

set(p,'LineWidth',2);
hold off

figure 
% policy after move
plotBoundary(NewPolicy{14},Qmax,Qmin,Xmax,Xmin,alpha,'fill')

hold on
p = plot([QVec(NB),QVec(NB)],[Xmin,Xmax],'k');

set(p,'LineWidth',2);
hold off

% buying profit
figure 

p = plot(XVec,BP);
titlestr = sprintf('Buying Profit with Initial Policy at q = %g',QVec(NB));
title(titlestr)
set(p,'LineWidth',3);
hold on
scatter(XVec(6),0,'MarkerEdgeColor','k','MarkerFaceColor','k');
scatter(XVec(3),BP(3),'MarkerEdgeColor','r','MarkerFaceColor','r')
scatter(XVec(8),BP(8),'MarkerEdgeColor','r','MarkerFaceColor','r')
hold off
ylabel('Buying Profit','fontsize',13)
xlabel('Log Price x','fontsize',13)

%% Decisions to make
% axis([Qmin,Qmax,Xmin,Xmax])
% xlabel('Volume in Storage q','Fontsize',13);
% ylabel('Log price x','Fontsize',13);
% title('Decisions to make','Fontsize',15)

%% 1 Dimension VS 2 Dimensions
