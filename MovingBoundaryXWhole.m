%% This function moves the boundary along X at the same time instead of doing it one by one.
%% MoveIndicator is odd means selling boundary should be moved, otherwise it is even.
function  [NewPolicy,Valuefunction,HoldingProfit,BuyingProfit,SellingProfit] = ...
    MovingBoundaryXWhole(M,U,alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,NumX,NumQ,InitialPolicy,BuyingType,BuyingCostPara,SellingType,SellingCostPara,MovingType,MoveIndicator)
%% Calculate BuyUpper, BuyLower and Sell from InitialPolicy
    BuyUpper = zeros(1,NumQ);
    BuyLower = zeros(1,NumQ);
    Sell = ones(1,NumQ)*(NumX + 1);
    
    for i = 1:NumQ
        if sum(InitialPolicy(i,:) == 1) ~= 0
            BuyUpper(i) = find(InitialPolicy(i,:) == 1,1,'last');
            BuyLower(i) = find(InitialPolicy(i,:) == 1,1,'first');
        end

        if sum(InitialPolicy(i,:) == 2) ~= 0
            Sell(i) = find(InitialPolicy(i,:) == 2,1,'first');
        end
    end
    
    
    
    
%% Calculate V
    % When there is no buying region
    % Maybe I need to separate the first iteration of movingboundary
    % and the else      
   
    if min(Sell) == NumX + 1 && min(BuyUpper) == 0
        error('Initial policy is holding all the time.');
    end

    [Valuefunction,HoldingProfit,BuyingProfit,SellingProfit] = ...
        SolvePDEWithTables(M,U,alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,NumX,NumQ,InitialPolicy,BuyingType,BuyingCostPara,SellingType,SellingCostPara);
%       [Valuefunction,HoldingProfit,BuyingProfit,SellingProfit] = ...
%         SolvePDENoDiscretization(M,U,alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,NumX,NumQ,InitialPolicy,BuyingType,BuyingCostPara,SellingType,SellingCostPara);
  

    
    %% Update the boundary and move selling boundary first
    if strcmp(MovingType,'SellingFirst')
        % If there is no buying region
        if max(BuyUpper) == 0        
            TotalSellingChange = 0;
            for i = 1 : NumQ
                % Move selling boundary first
                if SellingProfit(i,Sell(i)-1) > 0
                    Sell(i) = find(SellingProfit(i,:) == max(SellingProfit(i,:)),1,'last');
                    TotalSellingChange = TotalSellingChange + 1;
                end
            end

            % If we can not move any selling boundary move the buying one
            if TotalSellingChange == 0
                for i = 1: NumQ
                    % if there is no boundary at all move upper and lower to the
                    % same point
                    if BuyLower(i) == 0
                        if max(BuyingProfit(i,:)) > 0
                            BuyLower(i) = find(BuyingProfit(i,:) == max(BuyingProfit(i,:)),1,'first');
                            BuyUpper(i) = BuyLower(i);
                        end
                    else 
                        if max(BuyingProfit(i,1:BuyLower(i)-1)) > 0
                            BuyLower(i) = find(BuyingProfit(i,1:BuyLower(i)-1) == max(BuyingProfit(i,1:BuyLower(i)-1)),1,'first');
                        end

                        if max(BuyingProfit(i,BuyUpper(i)+1:Sell(i)-1)) > 0
                            BuyUpper(i) = find(BuyingProfit(i,BuyUpper(i)+1:Sell(i)-1) == max(BuyingProfit(i,BuyUpper(i)+1:Sell(i)-1)),1,'last') + BuyUpper(i);
                        end
                    end
                end
            end
        else
    
            % Move selling and buying at the same time
            if SellingProfit(i,Sell(i)-1) > 0
                Sell(i) = find(SellingProfit(i,:) == max(SellingProfit(i,:)),1,'last');
            end

            if max(BuyingProfit(i,1:BuyLower(i)-1)) > 0
                BuyLower(i) = find(BuyingProfit(i,1:BuyLower(i)-1) == max(BuyingProfit(i,1:BuyLower(i)-1)),1,'first');
            end

            if max(BuyingProfit(i,BuyUpper(i)+1:Sell(i)-1)) > 0
                BuyUpper(i) = find(BuyingProfit(i,BuyUpper(i)+1:Sell(i)-1) == max(BuyingProfit(i,BuyUpper(i)+1:Sell(i)-1)),1,'last') + BuyUpper(i);
            end
        end
    end
    
    %% Update the boundary alternatively
    
    if strcmp(MovingType,'Alternatively')
        %MoveIndicator is odd. Move selling boundary.
        if rem(MoveIndicator,2) ==  1
            for i = 1 : NumQ
                if SellingProfit(i,Sell(i)-1) > 0
                    Sell(i) = find(SellingProfit(i,:) == max(SellingProfit(i,:)),1,'last');
                end
            end
        %MoveIndicator is even. Move buying boundary.
        else
            for i = 1 : NumQ
                if BuyUpper(i) == 0
                    if max(BuyingProfit(i,:)) > 0
                        BuyLower(i) = find(BuyingProfit(i,:) == max(BuyingProfit(i,:)),1,'first');
                        BuyUpper(i) = BuyLower(i);
                    end
                else
                    if max(BuyingProfit(i,1:BuyLower(i)-1)) > 0
                        BuyLower(i) = find(BuyingProfit(i,1:BuyLower(i)-1) == max(BuyingProfit(i,1:BuyLower(i)-1)),1,'first');
                    end

                    if max(BuyingProfit(i,BuyUpper(i)+1:Sell(i)-1)) > 0
                        BuyUpper(i) = find(BuyingProfit(i,BuyUpper(i)+1:Sell(i)-1) == max(BuyingProfit(i,BuyUpper(i)+1:Sell(i)-1)),1,'last') + BuyUpper(i);
                    end
                end
            end
        end
    end
    
    %% Update the boundary together
    
    if strcmp(MovingType,'Together')   
        for i = 1 : NumQ
            if SellingProfit(i,Sell(i)-1) > 0
                Sell(i) = find(SellingProfit(i,:) == max(SellingProfit(i,:)),1,'last');
            end
            if BuyUpper(i) == 0
                if max(BuyingProfit(i,:)) > 0
                    BuyLower(i) = find(BuyingProfit(i,:) == max(BuyingProfit(i,:)),1,'first');
                    BuyUpper(i) = BuyLower(i);
                end
            else
                if max(BuyingProfit(i,1:BuyLower(i)-1)) > 0
                    BuyLower(i) = find(BuyingProfit(i,1:BuyLower(i)-1) == max(BuyingProfit(i,1:BuyLower(i)-1)),1,'first');
                end

                if max(BuyingProfit(i,BuyUpper(i)+1:Sell(i)-1)) > 0
                    BuyUpper(i) = find(BuyingProfit(i,BuyUpper(i)+1:Sell(i)-1) == max(BuyingProfit(i,BuyUpper(i)+1:Sell(i)-1)),1,'last') + BuyUpper(i);
                end
            end
        end
    end
    
    %% Generate new policy
    
    NewPolicy = zeros(NumQ,NumX);
    for i = 1 : NumQ
        if BuyLower(i) >0
            for j = 1 : BuyLower(i) -1 
                NewPolicy(i,j) = 0;
            end

            for j = BuyLower(i) : BuyUpper(i)
                NewPolicy(i,j) = 1;
            end
        end

        for j = BuyUpper(i)+1 : Sell(i)-1
            NewPolicy(i,j) = 0;
        end

        for j = Sell(i) : NumX
            NewPolicy(i,j) = 2;
        end
    end
    
%     dbstop in plotBoundary if 'MoveIndicator == 6'

end
       
    