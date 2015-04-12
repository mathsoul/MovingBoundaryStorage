% Function GenerateMCBuyOperator() generates the operator of buying for
% Markov Chain value/policy iteration.
% selling, and buying.
% Inputs: None
% Outputs:  
%           V(i,j,k) = operator_buy V + constant_buy

function [operator_buy,constant_buy] = GenerateMCBuyOperator()
    
    global Qmax Qmin Xmin Xmax Smin Smax NumX NumQ NumS BuyingType BuyingCostPara
    
    [operator_buy, constant_buy] = InitOperator();
    
    [indexMat, indexVecQ, indexVecX, indexVecS] = NodeIndex();
    
    dq = (Qmax - Qmin)/(NumQ - 1);
    dx = (Xmax - Xmin)/(NumX - 1);
    ds = (Smax - Smin)/(NumS - 1);
    
    QVec = Qmin:dq:Qmax;
    XVec = Xmin:dx:Xmax;
    SVec = Smin:ds:Smax;
    
    % V_q - ( e^(x + s) + lambda(q)) = 0
    % (V(i+1,j,k) - V(i,j,k))/dq = e^(x+s) + lambda(q)
    % V(i,j,k) = V(i+1,j,k) -(e^(x+s)dq + lambda(q)*dq)
    
    for ijk = 1:NumQ*NumX*NumS
        i = indexVecQ(ijk);
        j = indexVecX(ijk);
        k = indexVecS(ijk);
        
        q = QVec(i);
        x = XVec(j);
        s = SVec(k);
        
        if(~isOnUpperBorder(i,'Q'))
            operator_buy(ijk,indexMat(i+1,j,k)) = 1;
            constant_buy(ijk) = - (exp(x + s)*dq + AggregatedBuyingCost(BuyingType,q,q+dq,x+s,Qmax,Qmin,NumQ,BuyingCostPara));
        end
    end
end