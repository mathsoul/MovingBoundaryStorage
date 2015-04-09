% plotBoundary plots the boundary described in the buyingBoundary and sellingBoundary 
% function. 

function plotBoundary(Policy,Qmax,Qmin,Xmax,Xmin,alpha,type)
    [NumX,NumQ] = size(Policy);
    XVec = linspace(Xmin,Xmax,NumX);
    QVec = linspace(Qmin,Qmax,NumQ);
    
    BuyUBound = zeros(NumQ,1);
    BuyLBound = zeros(NumQ,1);
    
    SellUBound = (NumX+1)*ones(NumX,1);
    SellLBound = (NumX+1)*ones(NumX,1);
    
    for i = 1:NumQ
        if isempty(find(Policy(i,:) == 1,1,'first')) ~= 1
            BuyLBound(i) = find(Policy(i,:) == 1,1,'first');
        end
        if isempty(find(Policy(i,:) == 1,1,'last')) ~=1
            BuyUBound(i) = find(Policy(i,:) == 1,1,'last');
        end
        if isempty(find(Policy(i,:) == 2,1,'first')) ~= 1
            SellLBound(i) = find(Policy(i,:) == 2,1,'first');
        end
        if isempty(find(Policy(i,:) == 2,1,'last')) ~=1
            SellUBound(i) = find(Policy(i,:) == 2,1,'last');
        end
    end
    BuyMin = find(BuyUBound ~= 0,1,'first');
    BuyMax = find(BuyUBound ~= 0,1,'last');
    SellMin = find(SellLBound ~= NumX + 1,1,'first');
    
    BuyBound = [QVec(BuyMin:BuyMax),QVec(BuyMax:-1:BuyMin);...
        XVec(BuyUBound(BuyMin:BuyMax)),XVec(BuyLBound(BuyMax:-1:BuyMin))];
    SellBound = [QVec(NumQ:-1:SellMin),QVec(SellMin:NumQ);...
        XVec(SellUBound(NumQ:-1:SellMin)),XVec(SellLBound(SellMin:NumQ))];
    
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