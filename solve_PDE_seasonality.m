
function [VF,profit_hold,profit_sell,profit_buy] = solve_PDE_seasonality(Para)
Para = ParaVec();

[kappa,sigma,alpha,Qmax,Qmin,Xmax,Xmin,Smax,Smin,beta,NumX,NumQ,...
    NumS,BuyingType,BuyingCostPara,SellingType,SellingCostPara] = ...
    ParaDiv(Para);

%% Create the grid
dq = (Qmax - Qmin)/(NumQ - 1);       %Deltas
dx = (Xmax - Xmin)/(NumX - 1);
ds = (Smax - Smin)/(NumS - 1);

%% Create Node index

[G,Gi,Gj,Gk] = NodeIndex(NumQ,NumX,NumS);

%% Computation begins

[A,b,A_hold,A_sell,A_buy,b_hold,b_sell,b_buy]= operator(Para);

V = linsolve(A,b);

Valuefunction = reshape4disp(V,NumQ,NumX,NumS);




    



               