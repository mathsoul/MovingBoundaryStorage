% [VFPolicyIteration,PolicyPolicyIteration,NumIteration,HoldingProfit,BuyingProfit,SellingProfit] = 
% MarkovChain_with_seasonality(alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,NumX,NumQ,MaxIteration,
% ErrorTol,DiffLevel,BuyingType,BuyingCostPara,SellingType,SellingCostPara)%,type
% Computes the optimal value function via Policy Iteration. 
%
% The inputs are:
%   alpha is the long term mean value
%   beta is the discount factor
%   kappa is the rate the value reverts to the mean
%   sigma is the volatility 
%   Xmax is the maximum value of log price
%   Xmin is the minimum value of log price
%   Qmax is the maximum value of volume
%   Qmin is the minimum value of volume
%   NumX is the number of pieces used to discretize log price
%   NumQ is the number of pieces used to discretize volume
%   MaxIteration is the maximum number of iterations
%   ErrorTol is the tollerance of the error.
%
% and returns:
%   VFPI is the value function from policy iteration
%   PolicyPolicyIteration is the optimal policy of policy iteration
%   NumIteration is the number of iterations
%   HoldingProfit
%   BuyingProfit
%   SellingProfit


function [VF_PI,policy_PI,n_iter,profit_hold,profit_buy,profit_sell] = ...%,VFValueIteration,PolicyValueIteration
    MarkovChain_with_seasonality(alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,NumX,NumQ,MaxIteration,ErrorTol,DiffLevel,BuyingType,BuyingCostPara,SellingType,SellingCostPara)%,type
%% Create the grid

dq = (Qmax - Qmin)/(NumQ - 1);       %Deltas
dx = (Xmax - Xmin)/(NumX - 1);
ds = (Smax - Smin)/(NumS - 1);
XVec = Xmin:dx:Xmax;
QVec = Qmin:dq:Qmax;
SVec = Smin:ds:Smax;


%% Calculate the transaction cost

lambda = zeros(NumQ,NumX,NumS);
mu = zeros(NumQ,NumX,NumS);

for i = 1 : NumQ
    for j = 1 : NumX
        for k = 1:NumS
            lambda(i,j,k) = buyingCost(BuyingType,QVec(i),XVec(j)+SVec(k),Qmin,Qmax,BuyingCostPara);
            mu(i,j,k) = sellingCost(SellingType,QVec(i),XVec(j)+SVec(k),Qmin,Qmax,SellingCostPara);
        end
    end
end



%% Create Node index
[G,Gi,Gj,Gk] = NodeIndex(NumQ,NumX,NumS);

%% Computation begins



%Value iteration
MC{1}= zeros(NumQ*NumX,1);
run = 1; %converge or not
k = 1; %the number of iteration

while(run && k<MaxIteration)
%% Two variable MC{1} and MC{2}
    [MC{2},Policy] = max([Aholding * MC{1} + Cholding, Abuying * MC{1} + Cbuying, Aselling * MC{1} + Cselling],[],2);
    
    if(norm(MC{2}-MC{1}) < ErrorTol)
        run = 0;
    end
    MC{1} = MC{2}; 
    k = k+1;
end

VFValueIteration = MC{1};
PolicyValueIteration = Policy - 1;

clear Policy


%Policy iteration
Policy{1} = PolicyValueIteration;
run = 1; %converge or not
k = 1; %the number of iteration

while(run && k<MaxIteration)
    if(Policy{k}(ijk) == 1)
       A(ijk,:) = A_buy(ijk,:);
       b(ijk) = b_buy(ijk);
    elseif(Policy{k}(ijk) == 2)
        A(ijk,:) = A_sell(ijk,:);
        b(ijk) = b_sell(ijk);
    else
        A(ijk,:)=A_hold(ijk,:);
        b(ijk) = b_hold(ijk);
    end
    MC{k} = linsolve(eye(NumQ*NumX*NumS) - A,b); % I am not sure what this should be
    V_hold = A_hold * MC{k} + b_hold;
    V_buy = A_buy * MC{k} + b_buy;
    V_sell = A_sell * MC{k} + b_sell;
    [a,b] = max([V_hold,V_buy,V_sell],[],2);
    
    
    
    for ijk = 1:NumQ*NumX*NumS
        if(a(ijk)> ( 1+ DiffLevel)*MC{k}(ijk))
            Policy{k+1}(ijk,1) = b(ijk) - 1;
        else
            Policy{k+1}(ijk,1) = Policy{k}(ijk,1);
        end
    end

    if(norm(Policy{k+1}-Policy{k}) < ErrorTol)
        run = 0;
    end
    k = k+1;
end

VF_PI = MC{k-1};
policy_PI = Policy{k-1};
n_iter = k-1;

profit_hold = reshape4disp(V_hold,NumQ,NumX,NumS);

profit_buy = reshape4disp(V_hold,NumQ,NumX,NumS);

profit_sell = reshape4disp(V_hold,NumQ,NumX,NumS);

end












