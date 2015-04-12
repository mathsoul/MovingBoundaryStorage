% Function GenerateEquation() generates the linear equations for holding,
% selling, and buying.
% Inputs: None
% Outputs:  A x = b is the default equation
%           A_hold: holding A matrix
%           A_sell: selling A matrix
%           A_buy: buying A matrix
%           b_hold: holding b vector
%           b_sell: selling b vector
%           b_buy: buying b vector


function [A_hold,A_sell,A_buy,b_hold,b_sell,b_buy]= GenerateOperator()
    
   global kappa sigma alpha Qmax Qmin Xmin Xmax Smin Smax beta NumX NumQ...
        NumS BuyingType BuyingCostPara SellingType SellingCostPara
    
    [A_hold, b_hold] = InitEquation();
    [A_buy, b_buy] = InitEquation();
    [A_sell,b_sell] = InitEquation();
    
    [indexMat, indexVecQ, indexVecX, indexVecS] = NodeIndex();
    
    % These constants are used in computing hypergeometric functions
%     c1=beta/(2*kappa);
%     c2=kappa/sigma^2;

    for ijk = 1:NumQ*NumX*NumS
        i = indexVecQ(ijk);
        j = indexVecX(ijk);
        k = indexVecS(ijk);
        
        q = (i-1)*dq + Qmin;
        x = (j-1)*dx + Xmin;
        s = (k-1)*ds + Smin;

        %% second order discretize
        %holding 
        %The coefficient of V_xx V_x and V 
        coef_Vxx = 0.5*sigma^2/dx^2;
        coef_Vx = 0.5*kappa*(alpha-x)/dx;
        coef_V = (-2*coef_Vxx - beta);
        coef_Vs = 0.5*sqrt(1-s^2)/ds;

        if(ismember(j,2:(NumX-1)) && ismember(k,2:(NumS-1)))
            A_hold(ijk,indexMat(i,j+1,k)) = (coef_Vxx+coef_Vx)/coef_V;
            A_hold(ijk,indexMat(i,j-1,k)) = (coef_Vxx-coef_Vx)/coef_V;
            A_hold(ijk,indexMat(i,j,k+1)) = coef_Vs/coef_V;
            A_hold(ijk,indexMat(i,j,k-1)) = -coef_Vs/coef_V;
        elseif(j == 1)
            A_hold(ijk,ijk) = 1;
            A_hold(ijk,indexMat(i,j+1,k))= - 0.99;
            %-mchgu(c1,1/2,c2*(XVec(j)-alpha)^2)/mchgu(c1,1/2,c2*(XVec(j+1)-alpha)^2);
            % Kumar said this value doesn't matter that much
        elseif(j == NumX)
            A_hold(ijk,ijk) = 1;
            A_hold(ijk,indexMat(i,j-1,k))= 0.99;
            %-mchgu(c1,1/2,c2*(XVec(j)-alpha)^2)/mchgu(c1,1/2,c2*(XVec(j-1)-alpha)^2);
        elseif(k == 1)
            A_hold(ijk,indexMat(i,j+1,k)) = (coef_Vxx+coef_Vx)/(coef_V + 2*coef_Vs);
            A_hold(ijk,indexMat(i,j-1,k)) = (coef_Vxx-coef_Vx)/(coef_V + 2*coef_Vs);
            A_hold(ijk,indexMat(i,j,k+1)) = coef_Vs/(coef_V + 2*coef_Vs);
        elseif(k == NumS)
            A_hold(ijk,indexMat(i,j+1,k)) = (coef_Vxx+coef_Vx)/(coef_V + 2*coef_Vs);
            A_hold(ijk,indexMat(i,j-1,k)) = (coef_Vxx-coef_Vx)/(coef_V + 2*coef_Vs);
            A_hold(ijk,indexMat(i,j,k-1)) = coef_Vs/(coef_V + 2*coef_Vs);
        end


        % buying operator
        if(i < NumQ)
            A_buy(ijk,indexMat(i+1,j,k))= -1;
            b_buy(ijk) = - (exp(x + s)*dq + AggregatedBuyingCost(BuyingType,q,q+dq,x+s,Qmax,Qmin,NumQ,BuyingCostPara));
        end


        % selling operator
        if(1<i)
            A_sell(ijk,indexMat(i-1,j,k))=-1;
            b_sell(ijk) = exp(x + s) * dq - AggregatedSellingCost(SellingType,q-dq,q,x+s,Qmax,Qmin,NumQ,SellingCostPara);
        end
    end
end