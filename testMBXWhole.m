kappa = 3.4;
sigma =0.59;
alpha = 0.803;

Qmax=100;
Qmin=0;

% mark1
% Xmax= 2;
% Xmin = -0.6;
% beta = 0.05;
% 
% BuyingType = 'linear';
% BuyingCostPara = [0.02,0.05,-0.1]; 
% 
% SellingType = 'linear';
% SellingCostPara = [0.02,0.05,-0.1]; 

% 

%mark2
Xmax= 2.2;
Xmin = -4;

beta = 0.5;

BuyingType = 'reciprocal';
BuyingCostPara = [0.2,0,-0.2];

SellingType = 'reciprocal';
SellingCostPara = [0.2,0,-0.2];

%mark3
% Xmax= 2.2;
% Xmin = -4;
% 
% beta = 0.5;
% 
% BuyingType = 'reciprocal';
% BuyingCostPara = [0.2,0,-0.2];
% 
% SellingType = 'linear';
% SellingCostPara = [0.02,0,-0.2];

% mark4
% Xmax= 2.2;
% Xmin = -4;
% beta = 0.5;
% 
% BuyingType = 'linear';
% BuyingCostPara = [0.01,0.00,-0.1]; 
% 
% SellingType = 'linear';
% SellingCostPara = [0.02,0.05,-0.1]; 


NumX = 11;
NumQ = 11;



clear NewPolicy BuyingProfit SellingProfit


    
InitialPolicy = [zeros(1,NumX);zeros(NumQ-1,NumX-1),2*ones(NumQ-1,1)];

[NewPolicy{1},BuyingProfit{1},SellingProfit{1}] = ...
    MovingBoundaryXWhole(alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,NumX,NumQ,InitialPolicy,BuyingType,BuyingCostPara,SellingType,SellingCostPara);
Counter = 1;

while norm(NewPolicy{Counter}-InitialPolicy)>0
    InitialPolicy = NewPolicy{Counter};
    Counter = Counter +1;
   [NewPolicy{Counter},BuyingProfit{Counter},SellingProfit{Counter}] = ...
    MovingBoundaryXWhole(alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,NumX,NumQ,InitialPolicy,BuyingType,BuyingCostPara,SellingType,SellingCostPara);    
end


IMovingBoundaryXWhole = NewPolicy{Counter};



%adding constrains of when it is full no buying and when it is empty no
%selling.

Policy = IMovingBoundaryXWhole;

figure

hold on 

for i = 1 : NumQ
    for j = 1:NumX
        if(Policy(i,j)== 1)
            plot(Qmin + (i-1)*(Qmax-Qmin)/(NumQ-1),Xmin + (j-1)*(Xmax-Xmin)/(NumX-1),'o')
        elseif(Policy(i,j) == 2)
            plot(Qmin + (i-1)*(Qmax-Qmin)/(NumQ-1),Xmin + (j-1)*(Xmax-Xmin)/(NumX-1),'*')
        end        
    end
end

plot([Qmin,Qmax],[alpha,alpha])

axis([Qmin,Qmax,Xmin,Xmax])
% title({['Optimal policy from policy iteration with beta =', num2str(beta),' difference level =', num2str(DiffLevel)];['NumX = ', num2str(NumX); 'NumQ = ', num2str(NumQ)]});

% plot(QVec,log(CostRatio.*Qmax./QVec),'-')
hold off 