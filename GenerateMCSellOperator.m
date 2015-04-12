% Function GenerateMCSellOperator() generates the operator of selling for
% Markov Chain value/policy iteration.
% selling, and buying.
% Inputs: None
% Outputs:  
%           V(i,j,k) = operator_sell V + constant_sell

function [operator_sell,constant_sell] = GenerateMCSellOperator()
    
    global Qmax Qmin Xmin Xmax Smin Smax NumX NumQ NumS SellingType SellingCostPara
    
    [operator_sell, constant_sell] = InitOperator();
    
    [indexMat, indexVecQ, indexVecX, indexVecS] = NodeIndex();
    
    dq = (Qmax - Qmin)/(NumQ - 1);
    dx = (Xmax - Xmin)/(NumX - 1);
    ds = (Smax - Smin)/(NumS - 1);
    
    QVec = Qmin:dq:Qmax;
    XVec = Xmin:dx:Xmax;
    SVec = Smin:ds:Smax;
    % V_q - ( e^(x + s) - mu(q)) = 0
    % (V(i,j,k) - V(i-1,j,k))/dq = e^(x+s) - mu(q)
    % V(i,j,k) = V(i-1,j,k) + e^(x+s)dq - mu(q)*dq 
    
    for ijk = 1:NumQ*NumX*NumS
        i = indexVecQ(ijk);
        j = indexVecX(ijk);
        k = indexVecS(ijk);
        
        q = QVec(i);
        x = XVec(j);
        s = SVec(k);
        
        if(~isOnLowerBorder(i,'Q'))
            operator_sell(ijk,indexMat(i-1,j,k)) = 1;
            constant_sell(ijk) = exp(x + s) * dq - AggregatedSellingCost(SellingType,q-dq,q,x+s,Qmax,Qmin,NumQ,SellingCostPara);
        end
    end
end