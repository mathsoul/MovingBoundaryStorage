tic

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
% BuyingCostPara = [0.02,0.005,-0.1]; 
% 
% SellingType = 'linear';
% SellingCostPara = [0.02,0.005,-0.1]; 

% BuyingType = 'linear';
% BuyingCostPara = [0.1,0,-0.1];
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


NumX = 51;
NumQ = 51;

%% Move the boundary along X one by one.

% XVec = linspace(Xmin,Xmax,NumX);
% QVec = linspace(Qmin,Qmax,NumQ);
% 
% Counter = ones(1,NumQ);
% 
% clear NewPolicy BuyingProfit SellingProfit
% 
% for i = NumQ : -1 : 1
%     q = QVec(i);
%     
%     InitialPolicy = [zeros(1,NumX-1),2];
% 
%     [NewPolicy{i,1},BuyingProfit{i,1},SellingProfit{i,1}] = ...
%         MovingBoundaryAlongX(q,alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,NumX,NumQ,InitialPolicy,BuyingType,BuyingCostPara,SellingType,SellingCostPara);
% 
%     while norm(NewPolicy{i,Counter(i)}-InitialPolicy)>0
%         InitialPolicy = NewPolicy{i,Counter(i)};
%         Counter(i) = Counter(i) +1;
%        [NewPolicy{i,Counter(i)},BuyingProfit{i,Counter(i)},SellingProfit{i,Counter(i)}] = ...
%         MovingBoundaryAlongX(q,alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,NumX,NumQ,InitialPolicy,BuyingType,BuyingCostPara,SellingType,SellingCostPara);    
%     end
% end
% 
% IMovingBoundaryX = zeros(NumQ,NumX);
% 
% for i = 1:NumQ 
%     IMovingBoundaryX(i,:) = NewPolicy{i,Counter(i)};
% end
% 
% 
% %adding constrains of when it is full no buying and when it is empty no
% %selling.
% IMovingBoundaryX(NumQ,IMovingBoundaryX(NumQ,:)==1) =  0;
% IMovingBoundaryX(1,IMovingBoundaryX(1,:) == 2) =0;
% 
% Policy = IMovingBoundaryX;
% 
% figure
% 
% hold on 
% 
% for i = 1 : NumQ
%     for j = 1:NumX
%         if(Policy(i,j)== 1)
%             plot(Qmin + (i-1)*(Qmax-Qmin)/(NumQ-1),Xmin + (j-1)*(Xmax-Xmin)/(NumX-1),'o')
%         elseif(Policy(i,j) == 2)
%             plot(Qmin + (i-1)*(Qmax-Qmin)/(NumQ-1),Xmin + (j-1)*(Xmax-Xmin)/(NumX-1),'*')
%         end        
%     end
% end
% 
% plot([Qmin,Qmax],[alpha,alpha])
% 
% axis([Qmin,Qmax,Xmin,Xmax])
% % title({['Optimal policy from policy iteration with beta =', num2str(beta),' difference level =', num2str(DiffLevel)];['NumX = ', num2str(NumX); 'NumQ = ', num2str(NumQ)]});
% 
% % plot(QVec,log(CostRatio.*Qmax./QVec),'-')
% hold off
% 
% V4 = SolvePDENoDiscretization(alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,NumX,NumQ,IMovingBoundaryX,BuyingType,BuyingCostPara,SellingType,SellingCostPara);
% 
% [Aholding,Abuying,Aselling,Cholding,Cbuying,Cselling]= ...
%     OperatorGenerator(alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,NumX,...
%     NumQ,IMovingBoundaryX,BuyingType,BuyingCostPara,SellingType,SellingCostPara);
% 
%  V4cal = reshape((flipud(V4))',NumQ*NumX,1);
%     
%     Ih4 = flipud(reshape(Aholding*V4cal,NumQ,NumX)');
% 
%     Ib4 = flipud(reshape(-(Abuying*V4cal - Cbuying),NumQ,NumX)');
% 
%     Is4 = flipud(reshape(-(Aselling*V4cal - Cselling),NumQ,NumX)');
%     
% [V3,Ih3,Ib3,Is3] = SolvePDE(alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,NumX,NumQ,IMovingBoundaryX,BuyingType,BuyingCostPara,SellingType,SellingCostPara);
%  

%% Move the boundary along X at the same time
XVec = linspace(Xmin,Xmax,NumX);
QVec = linspace(Qmin,Qmax,NumQ);

Counter = 1;

clear NewPolicy BuyingProfit SellingProfit HoldingProfit Valuefunction m

M = zeros(1,NumX);
U = zeros(1,NumX);

a0 = beta/(2*kappa);
b0 = 1/2;

for i = 1:NumX
    M(i) = hypergeom(a0,b0,kappa/sigma^2*(XVec(i)-alpha)^2);
    U(i) = mchgu(a0,b0,kappa/sigma^2*(XVec(i)-alpha)^2);
end

InitialPolicy = [zeros(1,NumX);zeros(NumQ-1,NumX-1),2*ones(NumQ-1,1)];

[NewPolicy{1},Valuefunction{1},HoldingProfit{1},BuyingProfit{1},SellingProfit{1}] = ...
    MovingBoundaryXWhole(M,U,alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,NumX,NumQ,InitialPolicy,BuyingType,BuyingCostPara,SellingType,SellingCostPara,'Together');



while norm(NewPolicy{Counter}-InitialPolicy)>0
    InitialPolicy = NewPolicy{Counter};
    Counter = Counter +1;
   [NewPolicy{Counter},Valuefunction{Counter},HoldingProfit{Counter},BuyingProfit{Counter},SellingProfit{Counter}] = ...
    MovingBoundaryXWhole(M,U,alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,NumX,NumQ,InitialPolicy,BuyingType,BuyingCostPara,SellingType,SellingCostPara,'Together');    
end


IMovingBoundaryWholeX = NewPolicy{Counter};

toc

% HoldingProfitDisp = flipud(HoldingProfit{Counter}');
% IDisp = flipud(NewPolicy{Counter}');
% XVecDisp = fliplr(XVec);
% %adding constrains of when it is full no buying and when it is empty no
% %selling.
% 
% Policy = IMovingBoundaryWholeX;
% 
% figure
% 
% hold on 
% 
% for i = 1 : NumQ
%     for j = 1:NumX
%         if(Policy(i,j)== 1)
%             plot(Qmin + (i-1)*(Qmax-Qmin)/(NumQ-1),Xmin + (j-1)*(Xmax-Xmin)/(NumX-1),'o')
%         elseif(Policy(i,j) == 2)
%             plot(Qmin + (i-1)*(Qmax-Qmin)/(NumQ-1),Xmin + (j-1)*(Xmax-Xmin)/(NumX-1),'*')
%         end        
%     end
% end
% 
% plot([Qmin,Qmax],[alpha,alpha])
% 
% axis([Qmin,Qmax,Xmin,Xmax])
% % title({['Optimal policy from policy iteration with beta =', num2str(beta),' difference level =', num2str(DiffLevel)];['NumX = ', num2str(NumX); 'NumQ = ', num2str(NumQ)]});
% 
% % plot(QVec,log(CostRatio.*Qmax./QVec),'-')
% hold off
% 
% %% Plot the trajectory of moving boundary 
% for k  = 1:Counter
%     Policy = NewPolicy{k};
% 
% figure
% 
% hold on 
% 
% for i = 1 : NumQ
%     for j = 1:NumX
%         if(Policy(i,j)== 1)
%             plot(Qmin + (i-1)*(Qmax-Qmin)/(NumQ-1),Xmin + (j-1)*(Xmax-Xmin)/(NumX-1),'o')
%         elseif(Policy(i,j) == 2)
%             plot(Qmin + (i-1)*(Qmax-Qmin)/(NumQ-1),Xmin + (j-1)*(Xmax-Xmin)/(NumX-1),'*')
%         end        
%     end
% end
% 
% plot([Qmin,Qmax],[alpha,alpha])
% 
% axis([Qmin,Qmax,Xmin,Xmax])
% % title({['Optimal policy from policy iteration with beta =', num2str(beta),' difference level =', num2str(DiffLevel)];['NumX = ', num2str(NumX); 'NumQ = ', num2str(NumQ)]});
% m(k) = getframe;
% % plot(QVec,log(CostRatio.*Qmax./QVec),'-')
% hold off
% 
% 
% 
% % pause
% end
% 
% 
% close all
% 
% % Plot the final result from moving boundary whole X
% Policy = NewPolicy{k};
% 
% figure
% 
% hold on 
% 
% for i = 1 : NumQ
%     for j = 1:NumX
%         if(Policy(i,j)== 1)
%             plot(Qmin + (i-1)*(Qmax-Qmin)/(NumQ-1),Xmin + (j-1)*(Xmax-Xmin)/(NumX-1),'o')
%         elseif(Policy(i,j) == 2)
%             plot(Qmin + (i-1)*(Qmax-Qmin)/(NumQ-1),Xmin + (j-1)*(Xmax-Xmin)/(NumX-1),'*')
%         end        
%     end
% end
% 
% plot([Qmin,Qmax],[alpha,alpha])
% 
% axis([Qmin,Qmax,Xmin,Xmax])
% % title({['Optimal policy from policy iteration with beta =', num2str(beta),' difference level =', num2str(DiffLevel)];['NumX = ', num2str(NumX); 'NumQ = ', num2str(NumQ)]});
% m(k) = getframe;
% % plot(QVec,log(CostRatio.*Qmax./QVec),'-')
% hold off
% 
% %% Markov Chain method to achieve the optimal solution.
% tic
% MaxIteration = 1000;
% ErrorTol = 1;
% 
% % beta = [0.001,0.05,0.1,0.2];
% K=1;
% DiffLevel = 0 + (K-1)*0.001;
% 
% [V1,P1,N,Ih1,Ib1,Is1] = MarkovChain(alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,NumX,NumQ,MaxIteration,ErrorTol,DiffLevel,BuyingType,BuyingCostPara,SellingType,SellingCostPara);%V2,P2,
% V1 =flipud(reshape(V1,NumQ,NumX)');
% %         if(N>1)
% %             disp(N)
% %         end
% IMarkov = reshape(P1,NumQ,NumX);
% % [VMarkov,HoldingProfitMarkov,BuyingProfitMarkov,SellingProfitMarkov] = SolvePDENoDiscretization(alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,NumX,NumQ,IMarkov,BuyingType,BuyingCostPara,SellingType,SellingCostPara);
% 
% 
% 
% Policy = IMarkov;
% 
% figure
% 
% hold on 
% 
% for i = 1 : NumQ
%     for j = 1:NumX
%         if(Policy(i,j)== 1)
%             plot(Qmin + (i-1)*(Qmax-Qmin)/(NumQ-1),Xmin + (j-1)*(Xmax-Xmin)/(NumX-1),'o')
%         elseif(Policy(i,j) == 2)
%             plot(Qmin + (i-1)*(Qmax-Qmin)/(NumQ-1),Xmin + (j-1)*(Xmax-Xmin)/(NumX-1),'*')
%         end        
%     end
% end
% 
% plot([Qmin,Qmax],[alpha,alpha])
% 
% axis([Qmin,Qmax,Xmin,Xmax])
% title('Markov Chain')
% 
% m(Counter+1)= getframe;
% % plot(QVec,log(CostRatio.*Qmax./QVec),'-')
% hold off
% 
% movie2avi(m,'C:\Users\msbcr563\Google Drive\Longcode\Meeting91\Mark1_21by21.avi','compression','none','fps',1)
% 
% 
% % DiffMax = max(Valuefunction{Counter}-VMarkov);
% % LargerOrNot = Valuefunction{Counter} > VMarkov;
% % 
% % disp(LargerOrNot);
% % disp(DiffMax);
% 
% %% Plot the dynamics of SellingProfit and BuyingProfit
% 
% % clear mSellingProfit
% % 
% % hold on
% % for k = 1 : Counter
% %     l = (NumQ-1)/2;
% %     plot(XVec, SellingProfit{k}(l,:))
% % 
% %     mSellingProfit(k) = getframe;
% % end
% % hold off
% % 
% % movie2avi(mSellingProfit,'C:\Users\msbcr563\Google Drive\Longcode\Meeting91\Mark3SellingProfit_I_41by41.avi','compression','none','fps',1)
% % 
% % 
% % close all
% toc
% 
% UpOrNot = zeros(1,length(Valuefunction)-1);
% 
% for i = 1:(length(Valuefunction)-1)
%     UpOrNot(i) = all(all(Valuefunction{i+1}>Valuefunction{i})-0.001);
% end
    