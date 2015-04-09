

function  [Valuefunction,HoldingProfit, BuyingProfit, SellingProfit,NewPolicy,Bbound,Sbound] = ...
    CarefullyMovingBoundary(alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,NumX,NumQ,InitialPolicy,BuyingType,BuyingCostPara,SellingType,SellingCostPara)

    %% Generate the operator of calculating profits
    [Aholding,Abuying,Aselling,~,Cbuying,Cselling]= ...
    OperatorGenerator(alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,NumX,...
    NumQ,InitialPolicy,BuyingType,BuyingCostPara,SellingType,SellingCostPara);

    %% Calculate the buying bounds via Q axis

    Bbound = BuyingBoundsQaxis(InitialPolicy);

    %% Calculate the selling bounds via Q axis

    Sbound = SellingBoundsQaxis(InitialPolicy);
    
    %% Using non-discretized method to calculate
    Valuefunction = SolvePDENoDiscretization(alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,NumX,NumQ,InitialPolicy,BuyingType,BuyingCostPara,SellingType,SellingCostPara);    
    
    %% In order to calculate holding, buying and selling profit, valuefunction is transformed into a vector
    ValuefunctionCal = reshape((flipud(Valuefunction))',NumQ*NumX,1);
    
    HoldingProfit = reshape(-Aholding*ValuefunctionCal,NumQ,NumX);

    BuyingProfit = reshape(-(Abuying*ValuefunctionCal - Cbuying),NumQ,NumX);
    
    BuyingProfit(NumQ,:) = -inf; %When it is full, we can't buy.

    SellingProfit = reshape(-(Aselling*ValuefunctionCal - Cselling),NumQ,NumX);
    
    SellingProfit(1,:) = -inf; %When it is empty, we can't sell.
    
    %% Update the boundary
    Low = find(Bbound == 0,1,'first');
    Upper = find(Sbound == NumQ+1,1,'last');
    
    BuyingUnchange = 0;
    SellingUnchange = NumX +1;
    
    for j = 1 : NumX
        if BuyingProfit(Bbound(j)+1,j)> 0
            if j < Low
                Bbound(j) = find(BuyingProfit(:,j) == max(BuyingProfit(:,j)),1,'first');
            elseif j == Low
                if BuyingUnchange == j-1
                    Bbound(j) = find(BuyingProfit(:,j) == max(BuyingProfit(:,j)),1,'first');
                end
            end
        else
            if j < Low
                BuyingUnchange = j;
            end
        end
    end
    
    for j = NumX : -1 :1
        if SellingProfit(Sbound(j)-1,j)>0
            if j > Upper
                Sbound(j) = find(SellingProfit(:,j) == max(SellingProfit(:,j)),1,'first');
            elseif j ==Upper
                if SellingUnchange == j+1
                    Sbound(j) = find(SellingProfit(:,j) == max(SellingProfit(:,j)),1,'first');
                end
            end
        else
            if j > Upper
                SellingUnchange = j;
            end
        end
    end
    
    
    NewPolicy = PolicyGenerator(Bbound,Sbound,NumX,NumQ);
    
end
       
    