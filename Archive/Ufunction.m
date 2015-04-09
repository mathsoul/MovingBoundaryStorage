% NumT = 20;
% 
% XVec = linspace(-20,alpha,NumT); 
% UValue = zeros(1,NumT);
% 
% a0 = beta/(2*kappa);
% b0 = 1/2;
% 
% for i = 1:NumT
%     UValue(i) = mchgu(a0,b0,kappa/sigma^2*(XVec(i)-alpha)^2);
% end
% 
% 
% 
% plot(XVec,UValue)%./exp(XVec)
% 
% hold on
% 
% plot(XVec,exp(XVec))
% 
% %% Test whether the one dimensional moving boundary method is the one that give the optimal stopping problem.
% XVec = linspace(Xmin,Xmax,NumX);
% QVec = linspace(Qmin,Qmax,NumQ);
% 
% HoldingProfitBuy = zeros(1,NumX);
% HoldingProfitSell = zeros(1,NumX);
% 
% lambda = zeros(NumQ,NumX);
% mu = zeros(NumQ,NumX);
%     
% for i = 1 : NumQ
%     for j = 1:NumX
%         lambda(i,j) = buyingCost(BuyingType,QVec(i),XVec(j),Qmin,Qmax,BuyingCostPara);
%         mu(i,j) = sellingCost(SellingType,QVec(i),XVec(j),Qmin,Qmax,SellingCostPara);
%     end
% end
% 
% 
% i = NumQ - 6;
% 
% for j = 1:NumX
%     HoldingProfitBuy(j) = 0.5*sigma^2*exp(XVec(j))+kappa * (alpha - XVec(j))*exp(XVec(j)) - beta * (exp(XVec(j)) + lambda(i,j));
%     HoldingProfitSell(j) = 0.5*sigma^2*exp(XVec(j))+kappa * (alpha - XVec(j))*exp(XVec(j)) - beta * (exp(XVec(j)) - mu(i,j));
% end
% 
% 
% plot(XVec,HoldingProfitBuy)
% 
% figure
% 
% plot(XVec,HoldingProfitSell)

%% Test whether moving from right to left idea is right or not.

j =7;

Qend = 40;

LVq = zeros(1,Qend);

for i = 1 : Qend
    disp((0.5*sigma^2 + kappa *(alpha - XVec(j)) - beta) * exp(XVec(j)) );
    disp(beta * AggregatedBuyingCost(BuyingType,QVec(i),QVec(Qend),XVec(j),Qmax,Qmin,NumQ,BuyingCostPara));
    LVq(i) = (0.5*sigma^2 + kappa *(alpha - XVec(j)) - beta) * exp(XVec(j)) * (QVec(Qend) - QVec(i)) ...
        - beta * AggregatedBuyingCost(BuyingType,QVec(i),QVec(Qend),XVec(j),Qmax,Qmin,NumQ,BuyingCostPara);
end


plot(QVec(1:Qend),LVq);

%% 
% 
% V5 = SolvePDENoDiscretization(alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,NumX,NumQ,IMovingBoundaryTest,BuyingType,BuyingCostPara,SellingType,SellingCostPara);
% 
% [Aholding,Abuying,Aselling,Cholding,Cbuying,Cselling]= ...
%     OperatorGenerator(alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,NumX,...
%     NumQ,IMovingBoundaryX,BuyingType,BuyingCostPara,SellingType,SellingCostPara);
% 
%  V5cal = reshape((flipud(V5))',NumQ*NumX,1);
%     
%     Ih5 = flipud(reshape(Aholding*V5cal,NumQ,NumX)');
% 
%     Ib5 = flipud(reshape(-(Abuying*V5cal - Cbuying),NumQ,NumX)');
% 
%     Is5 = flipud(reshape(-(Aselling*V5cal - Cselling),NumQ,NumX)');

%% Test whether Vq satisfies the equation LVq = 0

Vqright = zeros(NumX,NumQ);
Vqleft = zeros(NumX,NumQ);
Vqright(:,1:NumQ-1) = (V4(:,2:NumQ) - V4(:,1:(NumQ-1)))/((Qmax-Qmin)/(NumQ-1));
Vqleft(:,2:NumQ) = (V4(:,2:NumQ) - V4(:,1:(NumQ-1)))/((Qmax-Qmin)/(NumQ-1));
Vqrightcal = reshape(flipud(Vqright)',NumQ*NumX,1);
Vqleftcal = reshape(flipud(Vqleft)',NumQ*NumX,1);

IhVqright = flipud(reshape(Aholding*Vqrightcal,NumQ,NumX)');
IhVqleft = flipud(reshape(Aholding*Vqleftcal,NumQ,NumX)');

