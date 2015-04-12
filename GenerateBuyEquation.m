% Function GenerateBuyEquation() generates the linear equations for holding,
% selling, and buying.
% Inputs: None
% Outputs:  A x = b is the default equation
%           A_buy 
%           b_buy 

function [A_buy,b_buy]= GenerateBuyEquation()
    
   global Qmax Qmin Xmin Xmax Smin Smax NumX NumQ NumS BuyingType BuyingCostPara
    
    [A_buy, b_buy] = InitEquation();
    
    [indexMat, indexVecQ, indexVecX, indexVecS] = NodeIndex();
    
    dq = (Qmax - Qmin)/(NumQ - 1);
    dx = (Xmax - Xmin)/(NumX - 1);
    ds = (Smax - Smin)/(NumS - 1);
    
    QVec = Qmin:dq:Qmax;
    XVec = Xmin:dq:Xmax;
    SVec = Smin:dq:Smax;

    for ijk = 1:NumQ*NumX*NumS
        i = indexVecQ(ijk);
        j = indexVecX(ijk);
        k = indexVecS(ijk);
        
        q = QVec(i);
        x = XVec(j);
        s = SVec(k);
        
        if(~isOnUpperBorder(i,'Q'))
            A_buy(ijk,indexMat(i+1,j,k))= -1;
            b_buy(ijk) = - (exp(x + s)*dq + AggregatedBuyingCost(BuyingType,q,q+dq,x+s,Qmax,Qmin,NumQ,BuyingCostPara));
        end
    end
end

    