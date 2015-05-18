
function [Aholding,Abuying,Aselling,Cholding,Cbuying,Cselling]= ...
    OperatorGenerator(U,alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,NumX,...
    NumQ,RegionIndicater,BuyingType,BuyingCostPara,SellingType,SellingCostPara)

%% Create the grid

dq=(Qmax-Qmin)/(NumQ-1);       %Deltas
dx=(Xmax-Xmin)/(NumX-1);
XVec=Xmin:dx:Xmax;
QVec=Qmin:dq:Qmax;

%% Calculate the transaction cost

% lambda = zeros(NumQ,NumX);
% mu = zeros(NumQ,NumX);
% 
% for i = 1 : NumQ
%     for j = 1 : NumX
%         lambda(i,j) = buyingCost(BuyingType,QVec(i),XVec(j),Qmin,Qmax,BuyingCostPara);
%         mu(i,j) = sellingCost(SellingType,QVec(i),XVec(j),Qmin,Qmax,SellingCostPara);
%     end
% end


%% Create Node index

G=repmat(0:NumQ:NumQ*NumX-1,NumQ,1)+repmat(1:NumQ,NumX,1)';        %Node index mat
Gi=reshape(repmat(1:NumQ,1,NumX),NumQ*NumX,1);                   % q index
Gj=reshape(repmat(1:NumX,NumQ,1),NumQ*NumX,1);                    % x index

%% Computation begins

% The discretized PDE is A V = C 
% A = zeros(NumQ*NumX);
% C = zeros(NumQ*NumX,1);

% Aholding, Abuying and Aselling are all operator matrix if all points are
% holding, buying and selling respectively. Similar for Cholding,
% Cbuying and Cselling.
Aholding=eye(NumQ*NumX);        
Abuying=eye(NumQ*NumX);      
Aselling=eye(NumQ*NumX);      

Cholding=zeros(NumQ*NumX,1);
Cbuying=zeros(NumQ*NumX,1);
Cselling=zeros(NumQ*NumX,1);


% Important constant (find out what are they?)

c1=beta/(2*kappa);
c2=kappa/sigma^2;

for ij=1:NumQ*NumX
    
    i=Gi(ij);j=Gj(ij);          %%Get i and j index
    q=(i-1)*dq+Qmin;
    x=(j-1)*dx+Xmin;

%% second order discretize
    %holding Region
    %The coefficient of V_xx V_x and V 
    CVxx=0.5*sigma^2/dx^2;
    CVx=0.5*kappa*(alpha-x)/dx;
    CV=(-2*CVxx - beta);
    
    if(j>1 && j<NumX)
        Aholding(ij,ij) = CV;
        Aholding(ij,G(i,j+1))=(CVxx+CVx);
        Aholding(ij,G(i,j-1))=(CVxx-CVx);
    elseif(j==1)
        Aholding(ij,G(i,j+1))= -U(j)/U(j+1);
%                 Aholding(ij,G(i,j+1))= -mchgu(c1,1/2,c2*(XVec(j)-alpha)^2)/mchgu(c1,1/2,c2*(XVec(j+1)-alpha)^2);
%         Aholding(ij,G(i,j+1))= -Hermite(-2*c1,sqrt(c2)*abs(XVec(j)-alpha))/Hermite(-2*c1,sqrt(c2)*abs(XVec(j+1)-alpha));
    elseif(j==NumX)
        Aholding(ij,G(i,j-1))= -U(j)/U(j-1);
%                 Aholding(ij,G(i,j-1))= -mchgu(c1,1/2,c2*(XVec(j)-alpha)^2)/mchgu(c1,1/2,c2*(XVec(j-1)-alpha)^2);
%         Aholding(ij,G(i,j-1))= -Hermite(-2*c1,sqrt(c2)*abs(XVec(j)-alpha))/Hermite(-2*c1,sqrt(c2)*abs(XVec(j-1)-alpha));
    end
    
    

    %buy
    if(i<NumQ)
        Abuying(ij,G(i+1,j))= -1;
%         Cbuying(ij)=-(exp(x)+lambda(i))*dq;
        Cbuying(ij) = - (exp(x)*dq + AggregatedBuyingCost(BuyingType,q,q+dq,XVec(j),Qmax,Qmin,NumQ,BuyingCostPara));
    else
        Cbuying(ij) = - inf;
    end
    
    
%     sell
    if(1<i)
        Aselling(ij,G(i-1,j))=-1;
%         Cselling(ij)=(exp(x)-mu(i))*dq;
        Cselling(ij) = exp(x) * dq - AggregatedSellingCost(SellingType,q-dq,q,XVec(j),Qmax,Qmin,NumQ,SellingCostPara);
    else
        Cselling(ij) = -inf;
    end
end