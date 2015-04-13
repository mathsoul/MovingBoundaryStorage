% Function GenerateSellEquation() generates the linear equations for holding,
% selling, and buying.
% Inputs: None
% Outputs:  A x = b is the default equation
%           A_sell 
%           b_sell 

function [A_sell,b_sell]= GenerateSellEquation()
    
   global Qmax Qmin Xmin Xmax Smin Smax NumX NumQ NumS SellingType SellingCostPara
    
    [A_sell, b_sell] = InitEquation();
    
    [indexMat, indexVecQ, indexVecX, indexVecS] = NodeIndex();
    
    dq = (Qmax - Qmin)/(NumQ - 1);
    dx = (Xmax - Xmin)/(NumX - 1);
    ds = (Smax - Smin)/(NumS - 1);
    
    QVec = Qmin:dq:Qmax;
    XVec = Xmin:dx:Xmax;
    SVec = Smin:ds:Smax;
    
    % V_q - ( e^(x + s) - mu(q)) = 0
    % (V(i,j,k) - V(i-1,j,k))/dq = e^(x+s) + mu(q)
    
    for ijk = 1:NumQ*NumX*NumS
        i = indexVecQ(ijk);
        j = indexVecX(ijk);
        k = indexVecS(ijk);
        
        q = QVec(i);
        x = XVec(j);
        s = SVec(k);
        
        if(~isOnLowerBorder(i,'Q'))
            A_sell(ijk,indexMat(i-1,j,k)) = -1;
            b_sell(ijk) = exp(x + s) * dq - AggregatedSellingCost(SellingType,q-dq,q,x+s,Qmax,Qmin,NumQ,SellingCostPara);
        else
            b_sell(ijk) = - Inf;
        end
    end
end

    