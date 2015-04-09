% [Valuefunction, RegionIndicater, HoldingProfit, BuyingProfit,
% SellingProfit] =
% SolvePDE(kappa,sigma,alpha,beta,Xmax,Xmin,Qmax,Qmin,NumQ,NumX,lambda,mu)
% Computes the value function via solving PDE numerically by discretizing
% it. Here the boundary is fixed and analytical which is discribed in
% another two functions buyingBoundary and sellingBoundary.
%
% The inputs are:
%   alpha is the long term mean value beta is the discount factor kappa is
%   the rate the value reverts to the mean sigma is the volatility Xmax is
%   the maximum value of log price Xmin is the minimum value of log price
%   Qmax is the maximum value of volume Qmin is the minimum value of volume
%   NumX is the number of pieces used to discretize log price NumQ is the
%   number of pieces used to discretize volume RegionIndicater is a NumX by
%   NumQ matrix where each element tells us which regiton the related point
%   belongs to. By using 1 representing buying, 2 representing selling and
%   0 representing holding. What's more, RegionIndicater is a NumQ by NumX
%   matrix. The value of volume is increasing along each column while the
%   log price is increasing along each row.
%
% and returns:
%   Valuefunction is a NumX by NumQ matrix. Each component is the value of
%   begining with related log price and volume. HoldingProfit is a NumX by
%   NumQ matrix where each element tells us what the instant profit is if
%   holding at the related point. BuyingProfit is a NumX by NumQ matrix
%   where each element tells us what the instant profit is if buying at the
%   related point. SellingProfit is a NumX by NumQ matrix where each
%   element tells us what the instant profit is if selling at the related
%   point.
function [V,HoldingProfit,BuyingProfit,SellingProfit] = ...
    SolvePDENoDiscretization(alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,NumX,NumQ,RegionIndicater,BuyingType,BuyingCostPara,SellingType,SellingCostPara)
%% Create the grid

dq=(Qmax-Qmin)/(NumQ-1);       %Deltas
dx=(Xmax-Xmin)/(NumX-1);
XVec=Xmin:dx:Xmax;
QVec=Qmin:dq:Qmax;

%% Calculate the holding bounds via X axis

[u,l] = HoldingBoundsXaxis(RegionIndicater);

%% Calculate the buying bounds via Q axis

Bbound = BuyingBoundsQaxis(RegionIndicater);
NumQMostBuying = max(Bbound);
% NumXMostBuying is the smallest number of log price where the amount we
% buy is most among all buyings.
NumXMostBuying = find(Bbound == NumQMostBuying,1,'first');

%% Calculate the selling bounds via Q axis

Sbound = SellingBoundsQaxis(RegionIndicater);

%% Calculate the parameters for HyperGeometric functions

a0 = beta/(2*kappa);
b0 = 1/2;


%% Calculate the transaction cost

lambda = zeros(NumQ,NumX);
mu = zeros(NumQ,NumX);

for i = 1 : NumQ
    for j = 1 : NumX
        lambda(i,j) = buyingCost(BuyingType,QVec(i),XVec(j),Qmin,Qmax,BuyingCostPara);
        mu(i,j) = sellingCost(SellingType,QVec(i),XVec(j),Qmin,Qmax,SellingCostPara);
    end
end


%% Create Node index

% G=reshape(1:NumQ*NumX,NumQ,NumX);        %Node index mat
% Gi=reshape(repmat(1:NumQ,1,NumX),NumQ*NumX,1);                   % q
% index Gj=reshape(repmat(1:NumX,NumQ,1),NumQ*NumX,1);                    %
% x index

%% Computation begins

% 
A = zeros(NumQ*3,NumQ*3);
b = zeros(NumQ*3,1);


for i = 1:NumQ
    A(3*(i-1)+1,3*(i-1)+1) = 1;
    A(3*(i-1)+1,3*(i-1)+2) = 2*gamma(b0)/gamma(a0-b0+1);
    A(3*(i-1)+1,3*(i-1)+3) = -1;
    b(3*(i-1)+1) = 0;
    if u(i) == NumX +1 %without upper bound
        A(3*(i-1)+2,3*(i-1)+1) = 1;
        b(3*(i-1)+2) = 0;
    else
%         if kappa/sigma^2*(XVec(u(i))-alpha)^2 > 19 %if the upper bound is too far away from alpha
%             if i == Sbound(u(i))
%                 A(3*(i-1)+2,3*(i-1)+1) = 1;
%                 A(3*(i-1)+2,3*(i-1)+2) = 0;
%                 A(3*(i-1)+2,3*(i-2)+1) = -1;
%                 A(3*(i-1)+2,3*(i-2)+2) = 0;
%                 b(3*(i-1)+2) = 0;
%     %             b(3*(i-1)+2) =
%     %             AggregatedSellingCost(QVec(i-1),QVec(i),Qmax);
%             else
%                 A(3*(i-1)+2,3*(i-1)+1) = 1;
%                 A(3*(i-1)+2,3*(i-1)+2) = 0;
%                 A(3*(i-1)+2,3*(Sbound(u(i))-1)+1) = -1;
%                 A(3*(i-1)+2,3*(Sbound(u(i))-1)+2) = 0;
%                 b(3*(i-1)+2) = 0;
%             end
%         else
            if i == Sbound(u(i))% (i,u(i)) is on the boundary
                A(3*(i-1)+2,3*(i-1)+1) = hypergeom(a0,b0,kappa/sigma^2*(XVec(u(i))-alpha)^2);
                A(3*(i-1)+2,3*(i-1)+2) = mchgu(a0,b0,kappa/sigma^2*(XVec(u(i))-alpha)^2);
                A(3*(i-1)+2,3*(i-2)+1) = -hypergeom(a0,b0,kappa/sigma^2*(XVec(u(i))-alpha)^2);
                A(3*(i-1)+2,3*(i-2)+2) = -mchgu(a0,b0,kappa/sigma^2*(XVec(u(i))-alpha)^2);
                b(3*(i-1)+2) = exp(XVec(u(i)))*dq - AggregatedSellingCost(SellingType,QVec(i-1),QVec(i),XVec(u(i)),Qmax,Qmin,NumQ,SellingCostPara);
            else
                A(3*(i-1)+2,3*(i-1)+1) = hypergeom(a0,b0,kappa/sigma^2*(XVec(u(i))-alpha)^2);
                A(3*(i-1)+2,3*(i-1)+2) = mchgu(a0,b0,kappa/sigma^2*(XVec(u(i))-alpha)^2);
                A(3*(i-1)+2,3*(Sbound(u(i))-1)+1) = -hypergeom(a0,b0,kappa/sigma^2*(XVec(u(i))-alpha)^2);
                A(3*(i-1)+2,3*(Sbound(u(i))-1)+2) = -mchgu(a0,b0,kappa/sigma^2*(XVec(u(i))-alpha)^2);
                b(3*(i-1)+2) = exp(XVec(u(i)))*(QVec(i)-QVec(Sbound(u(i)))) - AggregatedSellingCost(SellingType,QVec(Sbound(u(i))),QVec(i),XVec(u(i)),Qmax,Qmin,NumQ,SellingCostPara);
            end
%         end
    end
    if l(i) == 0 %without lower bound
        A(3*(i-1)+3,3*(i-1)+3) = 1;
        b(3*(i-1)+3) = 0;
    else
        if kappa/sigma^2*(XVec(l(i))-alpha)^2 > 30
             if i == Bbound(l(i))
                A(3*(i-1)+3,3*(i-1)+3) = 1;
                A(3*(i-1)+3,3*(i-1)+2) = 0;
                A(3*(i-1)+3,3*(i)+3) = -1;
                A(3*(i-1)+3,3*(i)+2) = 0;
                b(3*(i-1)+3) = 0;
            else
                A(3*(i-1)+3,3*(i-1)+3) = 1;
                A(3*(i-1)+3,3*(i-1)+2) = 0;
                A(3*(i-1)+3,3*(Bbound(l(i))-1)+3) = -1;
                A(3*(i-1)+3,3*(Bbound(l(i))-1)+2) = 0;
                b(3*(i-1)+3) = 0;
            end
        else
            if i == Bbound(l(i))
                A(3*(i-1)+3,3*(i-1)+3) = -hypergeom(a0,b0,kappa/sigma^2*(XVec(l(i))-alpha)^2);
                A(3*(i-1)+3,3*(i-1)+2) = mchgu(a0,b0,kappa/sigma^2*(XVec(l(i))-alpha)^2);
                A(3*(i-1)+3,3*(i)+3) = hypergeom(a0,b0,kappa/sigma^2*(XVec(l(i))-alpha)^2);
                A(3*(i-1)+3,3*(i)+2) = -mchgu(a0,b0,kappa/sigma^2*(XVec(l(i))-alpha)^2);
                b(3*(i-1)+3) = exp(XVec(l(i)))*dq + AggregatedBuyingCost(BuyingType,QVec(i),QVec(i+1),XVec(l(i)),Qmax,Qmin,NumQ,BuyingCostPara);
            else
                A(3*(i-1)+3,3*(i-1)+3) = -hypergeom(a0,b0,kappa/sigma^2*(XVec(l(i))-alpha)^2);
                A(3*(i-1)+3,3*(i-1)+2) = mchgu(a0,b0,kappa/sigma^2*(XVec(l(i))-alpha)^2);
                A(3*(i-1)+3,3*(Bbound(l(i))-1)+3) = hypergeom(a0,b0,kappa/sigma^2*(XVec(l(i))-alpha)^2);
                A(3*(i-1)+3,3*(Bbound(l(i))-1)+2) = -mchgu(a0,b0,kappa/sigma^2*(XVec(l(i))-alpha)^2);
                b(3*(i-1)+3) = exp(XVec(l(i)))*(QVec(Bbound(l(i)))-QVec(i))+ AggregatedBuyingCost(BuyingType,QVec(i),QVec(Bbound(l(i))),XVec(l(i)),Qmax,Qmin,NumQ,BuyingCostPara);
            end
        end
    end

    
    
end

X = linsolve(A,b);

V = zeros(NumQ,NumX);



for i = 1 : NumQ
    for j = NumXMostBuying: NumX
        if RegionIndicater(i,j) == 0
            if  XVec(j) > alpha 
                V(i,j) = X(3*(i-1)+1) * hypergeom(a0,b0,kappa/sigma^2*(XVec(j)-alpha)^2) + X(3*(i-1) +2)*mchgu(a0,b0,kappa/sigma^2*(XVec(j)-alpha)^2);
            else
                V(i,j) = X(3*(i-1)+3) * hypergeom(a0,b0,kappa/sigma^2*(XVec(j)-alpha)^2) - X(3*(i-1) +2)*mchgu(a0,b0,kappa/sigma^2*(XVec(j)-alpha)^2);
            end
%         elseif j == u(i)
%             V(i,j) = X(3*(i-1)+1) *
%             hypergeom(a0,b0,kappa/sigma^2*(XVec(j)-alpha)^2) + X(3*(i-1)
%             +2)*mchgu(a0,b0,kappa/sigma^2*(XVec(j)-alpha)^2);
%         elseif j == l(i)
%             V(i,j) = X(3*(i-1)+3) *
%             hypergeom(a0,b0,kappa/sigma^2*(XVec(j)-alpha)^2) - X(3*(i-1)
%             +2)*mchgu(a0,b0,kappa/sigma^2*(XVec(j)-alpha)^2);
        end
    end
end

if NumXMostBuying ~= 1
    for i = (NumQMostBuying + 1) : NumQ
        for j = 1: (NumXMostBuying -1)
            if RegionIndicater(i,j) == 0
                if  XVec(j) > alpha 
                    V(i,j) = X(3*(i-1)+1) * hypergeom(a0,b0,kappa/sigma^2*(XVec(j)-alpha)^2) + X(3*(i-1) +2)*mchgu(a0,b0,kappa/sigma^2*(XVec(j)-alpha)^2);
                else
                    V(i,j) = X(3*(i-1)+3) * hypergeom(a0,b0,kappa/sigma^2*(XVec(j)-alpha)^2) - X(3*(i-1) +2)*mchgu(a0,b0,kappa/sigma^2*(XVec(j)-alpha)^2);
                end
    %         elseif j == u(i)
    %             V(i,j) = X(3*(i-1)+1) *
    %             hypergeom(a0,b0,kappa/sigma^2*(XVec(j)-alpha)^2) +
    %             X(3*(i-1)
    %             +2)*mchgu(a0,b0,kappa/sigma^2*(XVec(j)-alpha)^2);
    %         elseif j == l(i)
    %             V(i,j) = X(3*(i-1)+3) *
    %             hypergeom(a0,b0,kappa/sigma^2*(XVec(j)-alpha)^2) -
    %             X(3*(i-1)
    %             +2)*mchgu(a0,b0,kappa/sigma^2*(XVec(j)-alpha)^2);
            end
        end
    end
end

%% Calculate the value function of buying region
%The part before curving back when the log price is decreasing
for j = NumXMostBuying :NumX
    for i = 1:NumQ
        if RegionIndicater(i,j) == 1
            if i == Bbound(j) %if it is on the boundary
                V(i,j) = V(i+1,j)-(exp(XVec(j))*dq + AggregatedBuyingCost(BuyingType,QVec(i),QVec(i+1),XVec(j),Qmax,Qmin,NumQ,BuyingCostPara));
            else
                V(i,j) = V(Bbound(j)+1,j) -(exp(XVec(j))*(QVec(Bbound(j)+1)-QVec(i)) + AggregatedBuyingCost(BuyingType,QVec(i),QVec(Bbound(j)+1),XVec(j),Qmax,Qmin,NumQ,BuyingCostPara));                
            end
        end
    end
end

%The curving back part
if NumXMostBuying ~=1
    XLimit = NumXMostBuying - 1; % Xlimit is the limit of log price we need to consider in this region
    for j = XLimit : -1 : 1
        V(NumQMostBuying,j) = mchgu(a0,b0,kappa/sigma^2*(XVec(j)-alpha)^2) / mchgu(a0,b0,kappa/sigma^2*(XVec(NumXMostBuying)-alpha)^2) * V(NumQMostBuying,NumXMostBuying);
    end
    for i = NumQMostBuying - 1 : -1 : max(Bbound(1),1)
        while XLimit > 0 && RegionIndicater(i,XLimit) == 1
            for k = 1:i
                V(k,XLimit) = V(i+1,XLimit) -(exp(XVec(XLimit))*(QVec(i+1)-QVec(k)) + AggregatedBuyingCost(BuyingType,QVec(k),QVec(i+1),XVec(j),Qmax,Qmin,NumQ,BuyingCostPara));
            end
            XLimit = XLimit - 1;
        end
        if XLimit>0 && RegionIndicater(i,XLimit) == 0
            for j = XLimit : -1 : 1
                V(i,j) = mchgu(a0,b0,kappa/sigma^2*(XVec(j)-alpha)^2) / mchgu(a0,b0,kappa/sigma^2*(XVec(XLimit +1)-alpha)^2) * V(i,XLimit+1);
            end
        end
    end
end


%% Calculate the value function of selling region
for j = 1 :NumX
    for i = 1:NumQ
        if RegionIndicater(i,j) == 2
            if i == Sbound(j) %if it is on the boundary
                V(i,j) = V(i-1,j)+(exp(XVec(j))*dq - AggregatedSellingCost(SellingType,QVec(i-1),QVec(i),XVec(j),Qmax,Qmin,NumQ,SellingCostPara));
            else
                V(i,j) = V(Sbound(j)-1,j) + (exp(XVec(j))*(QVec(i)- QVec(Sbound(j)-1)) - AggregatedSellingCost(SellingType,QVec(Sbound(j)-1),QVec(i),XVec(j),Qmax,Qmin,NumQ,SellingCostPara));                
            end
        end
    end
end

VCal = reshape(V,NumQ*NumX,1);
 
V = flipud(V');
% Valuefunction = flipud(reshape(V,NumQ,NumX)');

    [Aholding,Abuying,Aselling,~,Cbuying,Cselling]= ...
    OperatorGenerator(U,alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,NumX,...
    NumQ,RegionIndicater,BuyingType,BuyingCostPara,SellingType,SellingCostPara);



HoldingProfit = reshape(Aholding*VCal,NumQ,NumX);

BuyingProfit = reshape(-(Abuying*VCal - Cbuying),NumQ,NumX);

SellingProfit = reshape(-(Aselling*VCal - Cselling),NumQ,NumX);

BuyingProfit(NumQ,:) = -inf; %When it is full, we can't buy.
SellingProfit(1,:) = -inf; %When it is empty, we can't sell.
% end


% % Checks whether the dimension of lambda and mu are both NumQ function[]
% = checkDimensions(lambda, mu, NumQ)
% 
% if length(lambda) ~= NumQ
%     error('lambda does not have the same dimension as NumQ')
% end
% 
% if length(mu) ~= NumQ
%     error('mu does not have the same dimension as NumQ')
% end end











