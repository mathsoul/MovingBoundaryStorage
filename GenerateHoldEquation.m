% Function GenerateHoldEquation() generates the linear equations for holding,
% selling, and buying.
% Inputs: None
% Outputs:  A x = b is the default equation
%           A_hold 
%           b_hold 

function [A_hold,b_hold]= GenerateHoldEquation()
    
   global kappa sigma alpha Qmax Qmin Xmin Xmax Smin Smax beta NumX NumQ...
        NumS
    
    [A_hold, b_hold] = InitEquation();
    
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
        
        % For PDE a Vxx + b Vx + c V + d Vs = 0, we can discretize it as
        % following.
        
        % 1. Center discretization on both X and S dimension.
        
        % a (V(i,j+1,k) - 2V(i,j,k) + V(i,j-1,k))/(dx)^2 + 
        %   b(V(i,j+1,k) - V(i,j-1),k)/2dx +
        %   cV(i,j,k) + 
        %   d(V(i,j,k+1) - V(i,j,k-1))/2ds = 0
        
        % (a/(dx)^2 + b/2dx)) V(i,j+1,k) + (a/(dx)^2 - b/2dx) V(i,j-1,k) + 
        %   d/2ds V(i,j,k+1) - d/2ds V(i,j,k-1) = (2a/(dx)^2 - c) V(i,j,k)
        
        % 2. Center discretization on X and forward on S
        % (a/(dx)^2 + b/2dx)) V(i,j+1,k) + (a/(dx)^2 - b/2dx) V(i,j-1,k) + 
        %   d/ds V(i,j,k+1) = (2a/(dx)^2 - c + d/ds) V(i,j,k)
        
        % 3. Center discretization on X and backward on S
        % (a/(dx)^2 + b/2dx)) V(i,j+1,k) + (a/(dx)^2 - b/2dx) V(i,j-1,k) - 
        %   d/ds V(i,j,k-1) = (2a/(dx)^2 - c - d/ds) V(i,j,k)        
        
        % 1/2 sigma^2 Vxx + kappa (alpha - x) Vx - beta V + sqrt(1 - s^2) V_s = 0 
        
        coef_Vxx = 1/2*sigma^2;
        coef_Vx = kappa * (alpha - x);
        coef_V = -beta;
        coef_Vs = sqrt(1 - s^2);
        
        if(isInterior(j,'X') && isInterior(k,'S'))
            A_hold(ijk,indexMat(i,j+1,k)) = (coef_Vxx/dx^2 + coef_Vx/(2*dx))/(2*coef_Vxx/dx^2 - coef_V);
            A_hold(ijk,indexMat(i,j-1,k)) = (coef_Vxx/dx^2 - coef_Vx/(2*dx))/(2*coef_Vxx/dx^2 - coef_V);
            A_hold(ijk,indexMat(i,j,k+1)) = coef_Vs/(2*ds)/(2*coef_Vxx/dx^2 - coef_V);
            A_hold(ijk,indexMat(i,j,k-1)) = -coef_Vs/(2*ds)/(2*coef_Vxx/dx^2 - coef_V);
        elseif(isOnLowerBorader(j,'X'))
            A_hold(ijk,ijk) = 1;
            A_hold(ijk,indexMat(i,j+1,k))= - 0.99;
        elseif(isOnUpperBorader(j,'X'))
            A_hold(ijk,ijk) = 1;
            A_hold(ijk,indexMat(i,j-1,k))= 0.99;
        elseif(isOnLowerBorder(k,'S'))
            A_hold(ijk,indexMat(i,j+1,k)) = (coef_Vxx/dx^2 + coef_Vx/(2*dx))/(2*coef_Vxx/dx^2 - coef_V + coef_Vs/ds);
            A_hold(ijk,indexMat(i,j-1,k)) = (coef_Vxx/dx^2 - coef_Vx/(2*dx))/(2*coef_Vxx/dx^2 - coef_V + coef_Vs/ds);
            A_hold(ijk,indexMat(i,j,k+1)) = coef_Vs/ds/(2*coef_Vxx/dx^2 - coef_V + coef_Vs/ds);
        elseif(isOnUpperBorder(k,'S'))
            A_hold(ijk,indexMat(i,j+1,k)) = (coef_Vxx/dx^2 + coef_Vx/(2*dx))/(2*coef_Vxx/dx^2 - coef_V - coef_Vs/ds);
            A_hold(ijk,indexMat(i,j-1,k)) = (coef_Vxx/dx^2 - coef_Vx/(2*dx))/(2*coef_Vxx/dx^2 - coef_V - coef_Vs/ds);
            A_hold(ijk,indexMat(i,j,k-1)) = - coef_Vs/ds/(2*coef_Vxx/dx^2 - coef_V - coef_Vs/ds);
        end
    end
end