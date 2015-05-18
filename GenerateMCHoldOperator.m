% Function GenerateMCHoldOperator() generates the operator of holding for
% Markov Chain value/policy iteration.
% selling, and buying.
% Inputs: None
% Outputs:  
%           V(i,j,k) = operator_hold V + constant_hold

function [operator_hold,constant_hold]= GenerateMCHoldOperator()
    
   global kappa sigma alpha Qmax Qmin Xmin Xmax Smin Smax beta NumX NumQ...
        NumS
    
    [operator_hold, constant_hold] = InitOperator();
    
    [indexMat, indexVecQ, indexVecX, indexVecS] = NodeIndex();
    
    dq = (Qmax - Qmin)/(NumQ - 1);
    dx = (Xmax - Xmin)/(NumX - 1);
    ds = (Smax - Smin)/(NumS - 1);
    
    QVec = Qmin:dq:Qmax;
    XVec = Xmin:dx:Xmax;
    SVec = Smin:ds:Smax;

    for ijk = 1:NumQ*NumX*NumS
        i = indexVecQ(ijk);
        j = indexVecX(ijk);
        k = indexVecS(ijk);
        
        q = QVec(i);
        x = XVec(j);
        s = SVec(k);
                
        % 1/2 sigma^2 Vxx + kappa (alpha - x) Vx - beta V + V_s = 0 
        % All the elements in operator_hold should be non-negative
        
        % 1) If alpha > x
        %       (1/2 sigma^2/dx^2 + kappa(alpha - x)/dx) V(i,j+1,k) +
        %       (1/2 sigma^2/dx^2) V(i,j-1,k) + 
        %       (1/ds) V(i,j,k+1) =
        %       (beta + sigma^2/dx^2 + kappa(alpha-x)/dx + 1/ds) V(i,j,k)
        
        % 1) If alpha < x
        %       (1/2 sigma^2/dx^2) V(i,j+1,k) +
        %       (1/2 sigma^2/dx^2 + kappa(x - alpha)/dx) V(i,j-1,k) + 
        %       (1/ds) V(i,j,k+1) =
        %       (beta + sigma^2/dx^2 + kappa(x-alpha)/dx + 1/ds) V(i,j,k)
        
        c1 = 1/2*sigma^2/dx^2;
        c2 = kappa * abs(alpha - x)/dx;
        c3 = 1/ds;
        

        
        if(isInterior(j,'X'))
            if(~isOnUpperBorder(k,'S'))
                operator_hold(ijk,indexMat(i,j,k+1)) = c3;
            else
                operator_hold(ijk,indexMat(i,j,1)) = c3;
            end
            if( alpha > x)
                operator_hold(ijk,indexMat(i,j+1,k)) = (c1 + c2);
                operator_hold(ijk,indexMat(i,j-1,k)) = (c1);
            else
                operator_hold(ijk,indexMat(i,j+1,k)) = (c1);
                operator_hold(ijk,indexMat(i,j-1,k)) = (c1 + c2);
            end
            operator_hold(ijk,:) = operator_hold(ijk,:)/(beta + 2*c1 + c2 + c3);
        elseif(isOnLowerBorder(j,'X'))
            operator_hold(ijk,indexMat(i,j+1,k)) = 0.99;
        elseif(isOnUpperBorder(j,'X'))
            operator_hold(ijk,indexMat(i,j-1,k)) = 0.99;
        end
    end
end