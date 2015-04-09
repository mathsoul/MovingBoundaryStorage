

function  [NewPolicy,BuyingProfit,SellingProfit] = ...
    MovingBoundaryAlongX(q,alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,NumX,NumQ,InitialPolicy,BuyingType,BuyingCostPara,SellingType,SellingCostPara)%,MovingType
    %% Preparing
    XVec = linspace(Xmin,Xmax,NumX);
    a0 = beta/(2*kappa);
    b0 = 1/2;
    
    lambda = zeros(1,NumX);
    mu = zeros(1,NumX);
    
    
    for i = 1:NumX
        lambda(i) = buyingCost(BuyingType,q,XVec(i),Qmin,Qmax,BuyingCostPara);
        mu(i) = sellingCost(SellingType,q,XVec(i),Qmin,Qmax,SellingCostPara);
    end
    
    %% Calculate BuyUpper, BuyLower and Sell from InitialPolicy
    
    if sum(InitialPolicy == 1) ~= 0
        BuyUpper = find(InitialPolicy == 1,1,'last');
        BuyLower = find(InitialPolicy == 1,1,'first');
    else
        BuyUpper = 0;
        BuyLower = 0;
    end
    
    if sum(InitialPolicy == 2) ~= 0
        Sell = find(InitialPolicy == 2,1,'first');
    else
        Sell = NumX +1;
    end
    
    
    %% Calculate Vq
    Vq = zeros(1,NumX);
    
    % When there is no buying region
    while(BuyUpper == 0)
        if Sell == NumX + 1 
            error('Initial policy is holding all the time.');
        else
            xSell = XVec(Sell);
        end
    
        %% Vq in the selling region
        for i = Sell : NumX
            Vq(i) = exp(XVec(i)) - mu(i);
        end
        
        %% Vq in the holding region
        A = zeros(2,2);
        C = zeros(2,1);

        A(1,1) = hypergeom(a0,b0,kappa/sigma^2*(xSell-alpha)^2);
        A(1,2) = mchgu(a0,b0,kappa/sigma^2*(xSell-alpha)^2);
        A(2,1) = hypergeom(a0,b0,0);
        A(2,2) = 2*mchgu(a0,b0,0);

        C(1) = exp(xSell) - mu(Sell);
        C(2) = 0;

        Y = linsolve(A,C);

        for k = 1:Sell
            if XVec(k) > alpha
                Vq(k) = Y(1)*hypergeom(a0,b0,kappa/sigma^2*(XVec(k)-alpha)^2) + Y(2) * mchgu(a0,b0,kappa/sigma^2*(XVec(k)-alpha)^2);
            else
                Vq(k) = -Y(2) * mchgu(a0,b0,kappa/sigma^2*(XVec(k)-alpha)^2);
            end
        end

        %% Calculate selling profit and buying profit
        SellingProfit = -Vq + (exp(XVec)- mu);
        BuyingProfit = Vq - (exp(XVec) + lambda);
    
        %% Update the boundary
        % find the maximum and move
        if SellingProfit(Sell-1) > 0
            Sell = find(SellingProfit == max(SellingProfit),1,'last');
        else
            if max(BuyingProfit) > 0
                BuyLower = find(BuyingProfit == max(BuyingProfit),1,'first');
                BuyUpper = BuyLower;
            else
                break
            end
        end
        % find the positive and move(wrong)
%         if SellingProfit(Sell-1) > 0
%              Sell = find(SellingProfit > 0 ,1,'first');
%         else
%             if max(BuyingProfit) > 0
%                 BuyLower = find(BuyingProfit == max(BuyingProfit),1,'first');
%                 BuyUpper = BuyLower;
%             end
%         end
     end
        

        
  
    if BuyLower > 0
        xSell = XVec(Sell);
        xBuyUpper = XVec(BuyUpper);
        xBuyLower = XVec(BuyLower);


        %% Vq in the selling region
        for i = Sell : NumX
            Vq(i) = exp(XVec(i)) - mu(i);
        end

        %% Vq in the buying region
        for i = BuyLower : BuyUpper
            Vq(i) = exp(XVec(i)) + lambda(i);
        end

        %% Vq in the holding region which is above buying region

        A = zeros(3,3);
        C = zeros(3,1);

        A(1,1) = hypergeom(a0,b0,kappa/sigma^2*(xSell-alpha)^2);
        A(1,2) = mchgu(a0,b0,kappa/sigma^2*(xSell-alpha)^2);
        A(2,1) = hypergeom(a0,b0,0);
        A(2,2) = 2*mchgu(a0,b0,0);
        A(2,3) = -hypergeom(a0,b0,0);
        A(3,2) = - mchgu(a0,b0,kappa/sigma^2*(xBuyUpper-alpha)^2);
        A(3,3) = hypergeom(a0,b0,kappa/sigma^2*(xBuyUpper-alpha)^2);

        C(1) = exp(xSell) - mu(Sell);
        C(2) = 0;
        C(3) = exp(xBuyUpper)+lambda(BuyUpper);

        Y = linsolve(A,C);

        for k = BuyUpper+1:Sell-1
            if XVec(k) > alpha
                Vq(k) = Y(1)*hypergeom(a0,b0,kappa/sigma^2*(XVec(k)-alpha)^2) + Y(2) * mchgu(a0,b0,kappa/sigma^2*(XVec(k)-alpha)^2);
            else
                Vq(k) = -Y(2) * mchgu(a0,b0,kappa/sigma^2*(XVec(k)-alpha)^2) + Y(3)*hypergeom(a0,b0,kappa/sigma^2*(XVec(k)-alpha)^2);
            end
        end

        %% Vq in the holding region which is below buying region
        if XVec(BuyLower) > alpha
            error('Buying region is above alpha')
        else
            for k = 1 : BuyLower - 1
                Vq(k) = Vq(BuyLower)/mchgu(a0,b0,kappa/sigma^2*(xBuyLower-alpha)^2)*mchgu(a0,b0,kappa/sigma^2*(XVec(k)-alpha)^2);
            end
        end

        %% Calculate selling profit and buying profit
        SellingProfit = -Vq + (exp(XVec)- mu);
        BuyingProfit = Vq - (exp(XVec) + lambda);

        %% Update the boundary
        % find the maximum and move

        if SellingProfit(Sell-1) > 0
            Sell = find(SellingProfit == max(SellingProfit),1,'last');
        end

        if max(BuyingProfit(1:BuyLower-1)) > 0
            BuyLower = find(BuyingProfit(1:BuyLower-1) == max(BuyingProfit(1:BuyLower-1)),1,'first');
        end

        if max(BuyingProfit(BuyUpper+1:Sell-1)) > 0
            BuyUpper = find(BuyingProfit(BuyUpper+1:Sell-1) == max(BuyingProfit(BuyUpper+1:Sell-1)),1,'last') + BuyUpper;
        end

        %find the positive and move wrong!!!
    %     if SellingProfit(Sell-1) > 0
    %         Sell = find(SellingProfit >0 ,1,'first');
    %     end
    % 
    %     if max(BuyingProfit(1:BuyLower-1)) > 0
    %         BuyLower = find(BuyingProfit > 0,1,'first');
    %     end
    % 
    %     if max(BuyingProfit(BuyUpper+1:Sell-1)) > 0
    %         BuyUpper = find(BuyingProfit > 0,1,'last');
    %     end
    
    end
    
    %% Generate new policy
    NewPolicy = zeros(1,NumX);
    
    if BuyLower >0
        for i = 1 : BuyLower -1 
            NewPolicy(i) = 0;
        end

        for i = BuyLower : BuyUpper
            NewPolicy(i) = 1;
        end
    end

    for i = BuyUpper+1 : Sell-1
        NewPolicy(i) = 0;
    end

    for i = Sell : NumX
        NewPolicy(i) = 2;
    end
end
       
    