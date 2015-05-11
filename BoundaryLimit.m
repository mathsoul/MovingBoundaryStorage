% plotBoundary plots the boundary described in the buyingBoundary and sellingBoundary 
% function. 

function [buy_upper_limit, buy_lower_limit, sell_upper_limit, sell_lower_limit] = ...
    BoundaryLimit(policy)
    type = 'fill';
    global Qmax Qmin Xmax Xmin NumX NumQ alpha
    
    XVec = linspace(Xmin,Xmax,NumX);
    QVec = linspace(Qmin,Qmax,NumQ);
    
    buy_upper_limit = zeros(NumQ,1);
    buy_lower_limit = zeros(NumQ,1);
    
    sell_upper_limit = (NumX+1)*ones(NumX,1);
    sell_lower_limit = (NumX+1)*ones(NumX,1);
    
    for i = 1:NumQ
        if isempty(find(policy(i,:) == 1,1,'first')) ~= 1
            buy_lower_limit(i) = find(policy(i,:) == 1,1,'first');
        end
        if isempty(find(policy(i,:) == 1,1,'last')) ~=1
            buy_upper_limit(i) = find(policy(i,:) == 1,1,'last');
        end
        if isempty(find(policy(i,:) == 2,1,'first')) ~= 1
            sell_lower_limit(i) = find(policy(i,:) == 2,1,'first');
        end
        if isempty(find(policy(i,:) == 2,1,'last')) ~=1
            sell_upper_limit(i) = find(policy(i,:) == 2,1,'last');
        end
    end
    BuyMin = find(buy_upper_limit ~= 0,1,'first');
    BuyMax = find(buy_upper_limit ~= 0,1,'last');
    SellMin = find(sell_lower_limit ~= NumX + 1,1,'first');
    
    BuyBound = [QVec(BuyMin:BuyMax),QVec(BuyMax:-1:BuyMin);...
        XVec(buy_upper_limit(BuyMin:BuyMax)),XVec(buy_lower_limit(BuyMax:-1:BuyMin))];
    SellBound = [QVec(NumQ:-1:SellMin),QVec(SellMin:NumQ);...
        XVec(sell_upper_limit(NumQ:-1:SellMin)),XVec(sell_lower_limit(SellMin:NumQ))];
    
    if strcmp(type,'l')
        p = plot(SellBound(1,:),SellBound(2,:),'g',BuyBound(1,:),BuyBound(2,:),'r'...
            ,[Qmin,Qmax],[alpha,alpha],'b');
        set(p,'LineWidth',5);
    elseif strcmp(type,'fill')
        hold on
        p = plot(SellBound(1,:),SellBound(2,:),'g',BuyBound(1,:),BuyBound(2,:),'r'...
        ,[Qmin,Qmax],[alpha,alpha],'b');
        set(p,'LineWidth',5);
        fill(SellBound(1,:),SellBound(2,:),'g',BuyBound(1,:),BuyBound(2,:),'r');
        hold off
    end
    axis([Qmin,Qmax,Xmin,Xmax])
    xlabel('Volume in Storage q','Fontsize',13);
    ylabel('Log price x','Fontsize',13);
end