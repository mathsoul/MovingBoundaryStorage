%% Create the grid
kappa = 3.4;
sigma =0.59;
alpha = 0.803;

Qmax=100;
Qmin=0;

Xmax= 2.2;
Xmin = -4;

Smax = 1;
Smin = -1;

beta = 0.5;

NumX = 3;
NumQ = 3;
NumS = 3;

dq = (Qmax - Qmin)/(NumQ - 1);       %Deltas
dx = (Xmax - Xmin)/(NumX - 1);
ds = (Smax - Smin)/(NumS - 1);
XVec = Xmin:dx:Xmax;
QVec = Qmin:dq:Qmax;
SVec = Smin:ds:Smax;

BuyingType = 'reciprocal';
BuyingCostPara = [0.2,0,-0.2];

SellingType = 'reciprocal';
SellingCostPara = [0.2,0,-0.2];

%% Create Node index
NumQ = 3;
NumX = 3;
NumS = 3;

G = zeros(NumQ,NumX,NumS); %Node index mat
Gi = zeros(NumQ*NumX*NumS,1); % q index
Gj = zeros(NumQ*NumX*NumS,1); % x index
Gk = zeros(NumQ*NumX*NumS,1); % s index

for i = 1:NumQ
    for j = 1:NumX
        for k = 1:NumS
            G(i,j,k) = k + (j-1)*NumS + (i-1)*NumS*NumX; 
        end
    end
end

for l = 1:NumQ*NumX*NumS
    Gi(l) = ceil(l/(NumS * NumX));
    Gj(l) = ceil(rem(l,NumS*NumX)/NumS);
    Gk(l) = rem(l,NumS);
end

Gj(Gj == 0) = NumX;
Gk(Gk == 0) = NumS;

%% Computation begins

% The discretized PDE is A V = b 
A = zeros(NumQ*NumX*NumS);
b = zeros(NumQ*NumX*NumS,1);

A_hold = eye(NumQ*NumX*NumS);        
A_buy = eye(NumQ*NumX*NumS);      
A_sell = eye(NumQ*NumX*NumS);      

b_hold = zeros(NumQ*NumX*NumS,1);
b_buy = zeros(NumQ*NumX*NumS,1);
b_sell = zeros(NumQ*NumX*NumS,1);

% These constants are used in computing hypergeometric functions
c1=beta/(2*kappa);
c2=kappa/sigma^2;

for ijk = 1:NumQ*NumX*NumS
    i = Gi(ijk); j = Gj(ijk); k = Gk(ijk);
    q = (i-1)*dq + Qmin;
    x = (j-1)*dx + Xmin;
    s = (k-1)*ds + Smin;
    
    %% second order discretize
    %holding Region
    %The coefficient of V_xx V_x and V 
    coef_Vxx = 0.5*sigma^2/dx^2;
    coef_Vx = 0.5*kappa*(alpha-x)/dx;
    coef_V = (-2*coef_Vxx - beta);
    coef_Vs = 0.5*sqrt(1-s^2)/ds;
    
    if(ismember(j,2:(NumX-1)) && ismember(k,2:(NumS-1)))
        A_hold(ijk,G(i,j+1,k)) = (coef_Vxx+coef_Vx)/coef_V;
        A_hold(ijk,G(i,j-1,k)) = (coef_Vxx-coef_Vx)/coef_V;
        A_hold(ijk,G(i,j,k+1)) = coef_Vs/coef_V;
        A_hold(ijk,G(i,j,k-1)) = -coef_Vs/coef_V;
    elseif(j == 1)
        A_hold(ijk,ijk) = 1;
        A_hold(ijk,G(i,j+1,k))= - 0.99;
        %-mchgu(c1,1/2,c2*(XVec(j)-alpha)^2)/mchgu(c1,1/2,c2*(XVec(j+1)-alpha)^2);
        % Kumar said this value doesn't matter that much
    elseif(j == NumX)
        A_hold(ijk,ijk) = 1;
        A_hold(ijk,G(i,j-1,k))= 0.99;
        %-mchgu(c1,1/2,c2*(XVec(j)-alpha)^2)/mchgu(c1,1/2,c2*(XVec(j-1)-alpha)^2);
    elseif(k == 1)
        % Actually coef_Vs(i,j,k) is 0 when k = 1 or NumS
        A_hold(ijk,G(i,j+1,k)) = (coef_Vxx+coef_Vx)/(coef_V + 2*coef_Vs);
        A_hold(ijk,G(i,j-1,k)) = (coef_Vxx-coef_Vx)/(coef_V + 2*coef_Vs);
        A_hold(ijk,G(i,j,k+1)) = coef_Vs/(coef_V + 2*coef_Vs);
    elseif(k == NumS)
        A_hold(ijk,G(i,j+1,k)) = (coef_Vxx+coef_Vx)/(coef_V + 2*coef_Vs);
        A_hold(ijk,G(i,j-1,k)) = (coef_Vxx-coef_Vx)/(coef_V + 2*coef_Vs);
        A_hold(ijk,G(i,j,k-1)) = coef_Vs/(coef_V + 2*coef_Vs);
    end

    
    % buying region
    if(i < NumQ)
        A_buy(ijk,G(i+1,j,k))= -1;
        b_buy(ijk) = - (exp(x + s)*dq + AggregatedBuyingCost(BuyingType,q,q+dq,x,Qmax,Qmin,NumQ,BuyingCostPara));
    end
    
    
    % selling region
    if(1<i)
        A_sell(ijk,G(i-1,j,k))=-1;
        b_sell(ijk) = exp(x + s) * dq - AggregatedSellingCost(SellingType,q-dq,q,x,Qmax,Qmin,NumQ,SellingCostPara);
    end
    
    % assign region
    if(RegionIndicater(i,j,k) == 1)
       A(ijk,:) = A_buy(ijk,:);
       b(ijk) = b_buy(ijk);
    elseif(RegionIndicater(i,j,k) == 2)
        A(ijk,:) = A_sell(ijk,:);
        b(ijk) = b_sell(ijk);
    else
        A(ijk,:)=A_hold(ijk,:);
        b(ijk) = b_hold(ijk);
    end
end

V = linsolve(A,b);

Valuefunction = reshape(V,NumQ,NumX,NumS); %Test how to reshape it

HoldingProfit = reshape();



    



               