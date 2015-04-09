% [Valuefunction, RegionIndicater, HoldingProfit, BuyingProfit,
% SellingProfit] =
% SolvePDE(kappa,sigma,alpha,beta,Xmax,Xmin,Qmax,Qmin,NumQ,NumX,lambda,mu)
% Computes the value function via solving PDE numerically by discretizing
% it. Here the boundary is fixed and analytical which is discribed in
% another two functions buyingBoundary and sellingBoundary.
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
%   RegionIndicater is a NumX by NumQ matrix where each element tells us 
%   which regiton the related point belongs to. By using 1 representing 
%   buying, 2 representing selling and 0 representing holding.
%
% and returns:
%   Valuefunction is a NumX by NumQ matrix. Each component is the value of
%   begining with related log price and volume.
%   HoldingProfit is a NumX by NumQ matrix where each element tells us   
%   what the instant profit is if holding at the related point.
%   BuyingProfit is a NumX by NumQ matrix where each element tells us   
%   what the instant profit is if buying at the related point.
%   SellingProfit is a NumX by NumQ matrix where each element tells us   
%   what the instant profit is if selling at the related point.
function [Valuefunction,HoldingProfit, BuyingProfit, SellingProfit,A,Aholding,Abuying,Aselling,C,Cholding,Cbuying,Cselling] = ...
    SolvePDEForwardBackward(alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,NumX,NumQ,RegionIndicater,BuyingType,BuyingCostPara,SellingType,SellingCostPara)
%% Create the grid

dq=(Qmax-Qmin)/(NumQ-1);       %Deltas
dx=(Xmax-Xmin)/(NumX-1);
XVec=Xmin:dx:Xmax;
QVec=Qmin:dq:Qmax;

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

G=reshape(1:NumQ*NumX,NumQ,NumX);        %Node index mat
Gi=reshape(repmat(1:NumQ,1,NumX),NumQ*NumX,1);                   % q index
Gj=reshape(repmat(1:NumX,NumQ,1),NumQ*NumX,1);                    % x index

%% Computation begins

% The discretized PDE is A V = C 
A = zeros(NumQ*NumX);
C = zeros(NumQ*NumX,1);

% Aholding, Abuying and Aselling are all operator matrix if all points are
% holding, buying and selling respectively. Similar for Cholding,
% Cbuying and Cselling.
Aholding=-eye(NumQ*NumX);        
Abuying=-eye(NumQ*NumX);      
Aselling=-eye(NumQ*NumX);      

Cholding=zeros(NumQ*NumX,1);
Cbuying=zeros(NumQ*NumX,1);
Cselling=zeros(NumQ*NumX,1);


% Important constant (find out what are they?)

c1=beta/(2*kappa);
c2=kappa/sigma^2;

for ij=1:NumQ*NumX
    
    i=Gi(ij);j=Gj(ij);          %%Get i and j index
    q=(i-1)*dq+Qmin;
    x=(j-1)*dx+Xmin;
    
    %holding Region
    %The coefficient of V_xx V_x and V 
    CVxx=0.5*sigma^2/dx^2;
    CVx=kappa*(alpha-x)/dx;
    CV=(2*CVxx + abs(CVx)+beta);
    
    if(j>1 && j<NumX)
        if (x< alpha)
            Aholding(ij,G(i,j+1))=(CVxx+CVx)/CV;
            Aholding(ij,G(i,j-1))= CVxx/CV;
        else
            Aholding(ij,G(i,j+1))=CVxx/CV;
            Aholding(ij,G(i,j-1))=(CVxx-CVx)/CV;
        end
    elseif(j==1)
        Aholding(ij,G(i,j+1))= mchgu(c1,1/2,c2*(XVec(j)-alpha)^2)/mchgu(c1,1/2,c2*(XVec(j+1)-alpha)^2);
%         Aholding(ij,G(i,j+1))= Hermite(-2*c1,sqrt(c2)*abs(XVec(j)-alpha))/Hermite(-2*c1,sqrt(c2)*abs(XVec(j+1)-alpha));
%         Aholding(ij,G(i,j+1))= Hermite(-c1,sqrt(c2)*abs(XVec(j)-alpha))/Hermite(-2*c1,sqrt(c2)*abs(XVec(j+1)-alpha));
    elseif(j==NumX)
        Aholding(ij,G(i,j-1))= mchgu(c1,1/2,c2*(XVec(j)-alpha)^2)/mchgu(c1,1/2,c2*(XVec(j-1)-alpha)^2);
%         Aholding(ij,G(i,j-1))= Hermite(-2*c1,sqrt(c2)*abs(XVec(j)-alpha))/Hermite(-2*c1,sqrt(c2)*abs(XVec(j-1)-alpha));
%         Aholding(ij,G(i,j-1))= Hermite(-c1,sqrt(c2)*abs(XVec(j)-alpha))/Hermite(-2*c1,sqrt(c2)*abs(XVec(j-1)-alpha));
    end
    
    

    %buy
%     if(i<NumQ-1)
%         Abuying(ij,G(i+1,j))= 1;
% %         Cbuying(ij)=-(exp(x)+lambda(i))*dq;
%         Cbuying(ij) = - (exp(x)*dq + AggregatedBuyingCost(q,q+dq,Qmax,CostRatio));
%     elseif( i == NumQ-1)
%          Abuying(ij,G(i+1,j))= 1;
%          Cbuying(ij)=-(exp(x)+lambda(i))*dq;
%     end
%constant buying cost
    if(i<NumQ)
        Abuying(ij,G(i+1,j)) = 1;
%         Cbuying(ij) = -(exp(x)+lambda(i,j))*dq;
        Cbuying(ij) = (exp(x)*dq + AggregatedBuyingCost(BuyingType,QVec(i),QVec(i+1),XVec(j),Qmax,Qmin,NumQ,BuyingCostPara));
    end
    
    
%     sell
%     if(i>2)
%         Aselling(ij,G(i-1,j))=1;
% %         Cselling(ij)=(exp(x)-mu(i))*dq;
%         Cselling(ij) = exp(x) * dq - AggregatedSellingCost(q-dq,q,Qmax,CostRatio);
%     elseif(i ==2)
%         Aselling(ij,G(i-1,j))=1;
%         Cselling(ij)=(exp(x)-mu(i))*dq;
%     end
%constant selling cost
     if(i>1)
         Aselling(ij,G(i-1,j))= 1;
%          Cselling(ij) = (exp(x)-mu(i,j))*dq;
        Cselling(ij) = -(exp(x)*dq - AggregatedSellingCost(SellingType,QVec(i-1),QVec(i),XVec(j),Qmax,Qmin,NumQ,SellingCostPara));
     end
end


for ij=1:NumQ*NumX
    
    i=Gi(ij);j=Gj(ij);          %%Get i and j index
    
    if(RegionIndicater(i,j) == 1)
       A(ij,:) = Abuying(ij,:);
       C(ij) = Cbuying(ij);
    elseif(RegionIndicater(i,j) == 2)
        A(ij,:) = Aselling(ij,:);
        C(ij) = Cselling(ij);
    else
        A(ij,:)=Aholding(ij,:);
        C(ij) = Cholding(ij);
    end

end


V = linsolve(A,C);

Valuefunction = flipud(reshape(V,NumQ,NumX)');

HoldingProfit = reshape(Aholding*V,NumQ,NumX);

BuyingProfit = reshape((Abuying*V - Cbuying),NumQ,NumX);

SellingProfit = reshape((Aselling*V - Cselling),NumQ,NumX);

end


% % Checks whether the dimension of lambda and mu are both NumQ
% function[] = checkDimensions(lambda, mu, NumQ)
% 
% if length(lambda) ~= NumQ
%     error('lambda does not have the same dimension as NumQ')
% end
% 
% if length(mu) ~= NumQ
%     error('mu does not have the same dimension as NumQ')
% end
% end











