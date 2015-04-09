kappa = 3.4;
sigma =0.59;
alpha = 0.803;
beta = 0.05;

% syms x

Xmax= 1.3;
% Xmin=double(solve(1/2*sigma^2*exp(x)+kappa*(alpha-x-1)*exp(x) - beta*(exp(x))==0,x))-1;
% Xmin = 0.3;
Xmin = -0.2;

Qmax=100;
Qmin=0;

BuyingType = 'linear';
BuyingCostPara = [0.01,0.05,-0.1];

% BuyingType = 'linear';
% BuyingCostPara = [0.1,0,-0.1];
% 
% BuyingType = 'reciprocal';
% BuyingCostPara = [0.1,0,-0.1];

SellingType = 'linear';
SellingCostPara = [0.01,0.05,-0.1];

% SellingType = 'linear';
% SellingCostPara = [0.1,0,-0.1];
% 
% SellingType = 'reciprocal';
% SellingCostPara = [0.1,0,-0.1];
% 
% NumX = 41;
% NumQ = 41;

k = 1;

NumX = 21+(2^(k-1)-1)*20;
NumQ = 21+(2^(k-1)-1)*20;


clear Bbound Sbound HoldingProfit SellingProfit BuyingProfit Valuefunction


InitialPolicy = I;
% [Valuefunction{1},HoldingProfit,BuyingProfit{1},SellingProfit{1},NewPolicy,Bbound{1},Sbound{1}] = ...
%     CarefullyMovingBoundary(alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,NumX,NumQ,InitialPolicy,BuyingType,BuyingCostPara,SellingType,SellingCostPara);

[Valuefunction{1},HoldingProfit{1},BuyingProfit{1},SellingProfit{1},NewPolicy,Bbound{1},Sbound{1}] = ...
    MovingBoundary(alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,NumX,NumQ,InitialPolicy,BuyingType,BuyingCostPara,SellingType,SellingCostPara);



Counter = 1;

while norm(NewPolicy-InitialPolicy)>0
    Counter = Counter +1;
    InitialPolicy = NewPolicy;
%     [Valuefunction{Counter},HoldingProfit, BuyingProfit{Counter}, SellingProfit{Counter},NewPolicy,Bbound{Counter},Sbound{Counter}] = ...
%     CarefullyMovingBoundary(alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,NumX,NumQ,InitialPolicy,BuyingType,BuyingCostPara,SellingType,SellingCostPara);
    
    [Valuefunction{Counter},HoldingProfit{Counter}, BuyingProfit{Counter}, SellingProfit{Counter},NewPolicy,Bbound{Counter},Sbound{Counter}] = ...
    MovingBoundary(alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,NumX,NumQ,InitialPolicy,BuyingType,BuyingCostPara,SellingType,SellingCostPara);
    
end



result{1} = I;
result{2} = NewPolicy;
result{3} = Bbound;
result{4} = Sbound;
result{5} = HoldingProfit;
result{6} = SellingProfit;
result{7} = BuyingProfit;
result{8} = Valuefunction;

% XVec = linspace(Xmin,Xmax,NumX);
% 
% hold on
% for k = [1,19,20,length(result{3})]
% plot(XVec,result{8}{k}(NumX:-1:1,NumQ),'-');
% % plot(XVec,VMarkovChain(NumX:-1:1,NumQ),'--');
% % plot(XVec,VNoBuySellHighest(NumX:-1:1,NumQ),'-.');
% end
% hold off
% 
% title('q is 100');
% 

for k = [1:length(result{3})]%1,19,20,:length(result{3})
    hold on 
%     Policy = I;
      Policy = PolicyGenerator(result{3}{k},result{4}{k},NumX,NumQ);
    for i = 1 : NumQ
        for j = 1:NumX
            if(Policy(i,j)== 1)
                plot(Qmin + (i-1)*(Qmax-Qmin)/(NumQ-1),Xmin + (j-1)*(Xmax-Xmin)/(NumX-1),'o')
            elseif(Policy(i,j) == 2)
                plot(Qmin + (i-1)*(Qmax-Qmin)/(NumQ-1),Xmin + (j-1)*(Xmax-Xmin)/(NumX-1),'*')
            end        
        end
    end
%     
    plot([Qmin,Qmax],[alpha,alpha])

    axis([Qmin,Qmax,Xmin,Xmax])
    
    pause(0.2)
    
%     m(k) = getframe;
    hold off
end
% 
% % movie(m)
% 
% % XVec = linspace(Xmin,Xmax,NumX);
% % QVec = linspace(Qmin,Qmax,NumQ);
% % 
% % [MeshY,MeshX] = meshgrid(QVec,XVec);
% % mesh(MeshX,MeshY,result{5})
% % axis([0,100,0.3,1.3,-max(max(result{5})),max(max(result{5}))])
% 
% 
% 
% 
% 


%% Whether the valuefunction is increasing or not.

IncreaseRatio = zeros(length(result{8}) - 1,1);

for k = 1: (length(result{8}) - 1)
    IncreaseRatio(k) = sum(sum(result{8}{k+1}>result{8}{k}-0.001))/(NumX*NumQ);
end

plot( 1: length(result{8}) - 1,IncreaseRatio);


% plot(QVec,result{7}{1}(:,1))
% 
% hold on 
% 
% dq = (Qmax - Qmin)/(NumQ -1);
% 
% AggBuyingCost = zeros(1,NumQ-1);
% 
% for i = 1:NumQ-1
%     AggBuyingCost(i) = AggregatedBuyingCost(BuyingType,Qmin + (i-1)*dq,Qmin + i*dq,Xmin,Qmax,Qmin,NumQ,BuyingCostPara);
% end
% 
% plot(QVec(1:NumQ-1),-(exp(Xmin)*dq+AggBuyingCost),'--');
% 
% plot(QVec(1:NumQ-1),result{7}{1}(1:NumQ-1)+(exp(Xmin)*dq+AggBuyingCost),'-.');
% 
% 
% for i  = 2:Counter
%     disp(norm(result{3}{i} - result{3}{i-1},1))
% end
% 
% for i  = 2:Counter
%     disp(norm(result{4}{i} - result{4}{i-1},1))
% end