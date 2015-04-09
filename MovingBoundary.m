

function  [Valuefunction,HoldingProfit, BuyingProfit, SellingProfit,NewPolicy,Bbound,Sbound] = ...
    MovingBoundary(alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,NumX,NumQ,InitialPolicy,BuyingType,BuyingCostPara,SellingType,SellingCostPara)

    %% Generate the operator of calculating profits
%     [Aholding,Abuying,Aselling,~,Cbuying,Cselling]= ...
%     OperatorGenerator(alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,NumX,...
%     NumQ,InitialPolicy,BuyingType,BuyingCostPara,SellingType,SellingCostPara);

    %% Calculate the buying bounds via Q axis

    Bbound = BuyingBoundsQaxis(InitialPolicy);

    %% Calculate the selling bounds via Q axis

    Sbound = SellingBoundsQaxis(InitialPolicy);
    
%     %% Using non-discretized method to calculate
    Valuefunction = SolvePDENoDiscretization(alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,NumX,NumQ,InitialPolicy,BuyingType,BuyingCostPara,SellingType,SellingCostPara);    
    
    %% In order to calculate holding, buying and selling profit, valuefunction is transformed into a vector
    ValuefunctionCal = reshape((flipud(Valuefunction))',NumQ*NumX,1);
    
    HoldingProfit = reshape(-Aholding*ValuefunctionCal,NumQ,NumX);

    BuyingProfit = reshape(-(Abuying*ValuefunctionCal - Cbuying),NumQ,NumX);
    
    BuyingProfit(NumQ,:) = -inf; %When it is full, we can't buy.

    SellingProfit = reshape(-(Aselling*ValuefunctionCal - Cselling),NumQ,NumX);
    
    SellingProfit(1,:) = -inf; %When it is empty, we can't sell.

%% Using forward and backward method to calculate 

%     [Valuefunction,HoldingProfit,BuyingProfit,SellingProfit] = ...
%     SolvePDEForwardBackward(alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,NumX,NumQ,InitialPolicy,BuyingType,BuyingCostPara,SellingType,SellingCostPara);
% 
% BuyingProfit(NumQ,:) = -inf; %When it is full, we can't buy.
% SellingProfit(1,:) = -inf; %When it is empty, we can't sell.
    
    %% Update the boundary
    
    
    for j = 1 : NumX
        if BuyingProfit(Bbound(j)+1,j)> 0
            Bbound(j) = find(BuyingProfit(:,j) == max(BuyingProfit(:,j)),1,'first');
        end
        if SellingProfit(Sbound(j)-1,j)>0
            Sbound(j) = find(SellingProfit(:,j) == max(SellingProfit(:,j)),1,'first');
        end
    end
    
    NewPolicy = PolicyGenerator(Bbound,Sbound,NumX,NumQ);
    
end
       
    