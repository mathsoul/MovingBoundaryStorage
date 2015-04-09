

%set parallel options
% matlabpool(7)


tic
%% Parameters from Stathis' paper


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


% NumX = 41;
% NumQ = 41;


%% Parameters from Kumar's intial code


% sigma=.334;  %vol
% kappa=.301;  % rev speed
% gamma=3.093; %rev level for (ln S)
% alpha=(gamma - sigma^2/(2*kappa));  %risk neutral rev level for X
% beta=.3;
% 
% Xmax=3.6;
% Xmin=0.2; %log(Smin);
% 
% Qmin=0;
% Qmax=100;
% 
% NumX =21;
% NumQ =21;


%% 
%  plotBoundary(alpha,beta,kappa,sigma,Qmax,Qmin,Xmax,Xmin)
%  title('Boundary for the realistic parameters')
% 
% plot(logPriceGenerator(Xmax,alpha,kappa,sigma,tau,T/tau))
% 
% hold on 
% plotGrids(Xmax,Xmin,Qmax,Qmin,NumX,NumQ)
% hold off
% 
% 
% Num = [51,51,51];
% V = cell(3,3);

% Cond = zeros(3,3); %conditional number of a matrix

%% If analytical boundary is given in the sellingBoundary and buyingBoundary functions

% I = RegionIndicater(alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,NumX,NumQ); 
% 
% % for i = 1:3
% %     for j = 1:3
% %         NumX = Num(i);
% %         NumQ = Num(j);
% [V1,Ih,Ib,Is,A,Aholding,Abuying,Aselling,C,Cholding,Cbuying,Cselling] = SolvePDE(alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,NumX,NumQ,I);
%         V(i,j) = {V1};
% %         Cond(i,j) = cond(A);
%     end
% end


%% Simulation
% tau = 0.01;
% T = 100;
% N = 1000;
% % 
% V2 = zeros(NumX,NumQ);
% SD = zeros(NumX,NumQ);
% BT = zeros(NumX,NumQ);
% ST = zeros(NumX,NumQ);
% NS = zeros(NumX,NumQ);
% TB = zeros(NumX,NumQ);
% TS = zeros(NumX,NumQ);
% QInitial = zeros(NumX,NumQ);
% VInitial = zeros(NumX,NumQ);
% BInitial = zeros(NumX,NumQ);
% SInitial = zeros(NumX,NumQ);
% % 
% % 
% % 
% % % hold on
% % 
% % 
% % 
% % 
% % 
% for i = 1%[(1+NumQ)/2]%1,,1,NumQ
%     for j = 1:NumX%(NumX+1)/2,
%         x = Xmin + (j-1)*(Xmax-Xmin)/(NumX-1);
%         q = Qmin+ (i-1)*(Qmax-Qmin)/(NumQ-1);
%         [V2(NumX+1-j,i),SD(NumX+1-j,i),BT(NumX+1-j,i),ST(NumX+1-j,i),NS(NumX+1-j,i),TB(NumX+1-j,i),TS(NumX+1-j,i)]...
%            = ParSimulation(x,q,alpha,beta,kappa,sigma,Qmax,Qmin,T,N,tau);
% %        [QInitial(NumX+1-j,i),VInitial(NumX+1-j,i),~,~,...
% %         ~,BInitial(NumX+1-j,i),SInitial(NumX+1-j,i)] = ...
% %         InitialTransaction(x,q,alpha,beta,kappa,sigma,Qmax,Qmin,T,N,tau);
% %         
% %         VTime = ValueSimTime(x,q,alpha,beta,kappa,sigma,Qmax,Qmin,T,N,tau);
% % % 
% %         figure
% % 
% %         plot(1:T,VTime((1:T)/tau))
% %         title({['Valuefunction with respect to time'];['starting point is (' num2str(q) ',' num2str(x) ')']})
% % %         
% %         figure
% %         PlotOnePath(x,q,alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,T,tau,100)
% %         title({['One path'];['starting point is (' num2str(q) ',' num2str(x) ')'];[num2str(100) ' points ', num2str(1) ' unit time ']})
%     end
% end
% % Vdot1 = V2;
% % Vdot01 = V2;
% % Vdot001 =V2;
% % hold on
% % 
% % plot(linspace(Xmin,Xmax,NumX),flipud(Vdot1(:,NumQ)))
% % 
% % plot(linspace(Xmin,Xmax,NumX),flipud(Vdot01(:,NumQ)),'-.')
% % % 
% % plot(linspace(Xmin,Xmax,NumX),flipud(Vdot001(:,NumQ)),'--')
% % 
% % plot(linspace(Xmin,Xmax,NumX),flipud(V1(:,NumQ)),'*')
% % hold off
% % % 
% % % hold on
% % % 
% % % plot(linspace(Xmin,Xmax,NumX),flipud(V1(:,NumQ)))
% % % 
% % % scatter(linspace(Xmin,Xmax,NumX),flipud(SD(:,NumQ)))
% % % 
% % % title({['q =100, dashed line is from simulation, time limit is' num2str(T)]; ['Cirles are standard deviation'];['the length of each step is' num2str(tau)];['the number of simulation is'  num2str(N)]})
% % % 
% % % hold off
% % 
% % % hold off
% % 
% % V = reshape(flipud(V2)',NumX*NumQ,1);
% % 
% % Profit = flipud(reshape(-(A*V-C),NumQ,NumX)');
% % 
% % % HoldingProfit = flipud(reshape(-Aholding*V,NumQ,NumX)');
% % % 
% % % BuyingProfit = flipud(reshape(-(Abuying*V - Cbuying),NumQ,NumX)');
% % % 
% % % SellingProfit = flipud(reshape(-(Aselling*V - Cselling),NumQ,NumX)');
% % 
% % % ratio = V2(1,:)./V2(2,:);
% % % Diff = V1-V2;
% % 
% for n = [NumQ];%1,(NumQ+1)/2,
% 
%     plotValueFunction(n,V1,V2,SD,Xmax,Xmin,Qmax,Qmin,T,tau,N);
% % 
% % %     plotTransactionCosts(n,TB,TS,Xmax,Xmin,Qmax,Qmin,T,tau,N)
% % %     
% % %     plotTransactionTimes(n,BT,ST,Xmax,Xmin,Qmax,Qmin,T,tau,N)
% % %     
% % %     plot(linspace(Xmin,Xmax,NumX),fliplr(V2(:,n)+TB(:,n)+TS(:,n)));
% % %     
% %    
% %     
% % %     hold on 
% % % 
% % %     scatter(linspace(Xmin,Xmax,NumX),fliplr(Profit(:,n))')
% % % 
% % %     hold off
% % 
% end
% % 
% XVec = linspace(Xmin,Xmax,NumX);
% % 
% % % plotyy(XVec,flipud(V1(:,NumQ)),XVec,flipud(V2(:,NumQ)))
% % % 
% % % title('Value function from two methods with q = 100')
% % 
% % plot(XVec,flipud(V1(:,NumQ)),'-',XVec,flipud(Profit(:,NumQ)),'--')
% % 
% % title({['Value function from solving PDE and residuals with q = 100']; ['Dashed line is the residual']})
% % 
% % % IntialRatio = VInitial./V1 .*(I==2);
% % % 
% % % scatter(linspace(Xmin,Xmax,NumX),flipud(IntialRatio(:,NumQ)))
% % % 
% % % title({['Value after intial transaction over value of value function'],['q is ' num2str(Qmax)]})
% % % 
% % % figure
% % % 
% % % scatter(linspace(Xmin,Xmax,NumX),flipud(IntialRatio(:,(1+NumQ)/2)))
% % % 
% % % title({['Value after intial transaction over value of value function'],['q is ' num2str((Qmin+Qmax)/2)]})
% % 
% % %     X = linspace(Xmin,Xmax,NumX);
% % %     Q = linspace(Qmin,Qmax,NumQ);
% % %     
% % %     QGrid1 = ones(NumX,1)*Q;
% % %     XGrid1 = X'*ones(1,NumQ);
% % %     
% % % 
% % % scatter3(reshape(QGrid1,1,NumX*NumQ),reshape(XGrid1,1,NumX*NumQ),reshape(IntialRatio,1,NumX*NumQ));
% % 
% 
% TC = TB +TS;
% P = TC + V2;
% 
% plot(XVec,flipud(TC(:,1)),'-',XVec,flipud(P(:,1)),'--')
% 
% title({['Total transaction cost and profit from buying low and selling high. q = 0']; ['Dashed line is the profit']})

% toc


%% Markov Chain method to achieve the optimal solution.
MaxIteration = 1000;
ErrorTol = 1;

% beta = [0.001,0.05,0.1,0.2];
k=1;


DiffLevel = 0 + (k-1)*0.001;
        
for k = 2:5

    NumX = 21+(k-1)*20;
    NumQ = 21+(k-1)*20;
    XVec = linspace(Xmin,Xmax,NumX);
    QVec = linspace(Qmin,Qmax,NumQ);    
    [V1,P1,N,Ih1,Ib1,Is1] = MarkovChain(alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,NumX,NumQ,MaxIteration,ErrorTol,DiffLevel,BuyingType,BuyingCostPara,SellingType,SellingCostPara);%V2,P2,
    V1 =flipud(reshape(V1,NumQ,NumX)');
%         if(N>1)
%             disp(N)
%         end
    IMarkov = reshape(P1,NumQ,NumX);
%     HoldingMinusBuying = Ih1-Ib1;
%     [V2,Ih,Ib,Is,A,Aholding,Abuying,Aselling,C,Cholding,Cbuying,Cselling] = SolvePDE(alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,NumX,NumQ,IMarkov,BuyingType,BuyingCostPara,SellingType,SellingCostPara);
% % % 
%     [Aholding,Abuying,Aselling,Cholding,Cbuying,Cselling]= ...
%     OperatorGenerator(alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,NumX,...
%     NumQ,IMarkov,BuyingType,BuyingCostPara,SellingType,SellingCostPara);
% 
% 
%     V3 = SolvePDENoDiscretization(alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,NumX,NumQ,IMarkov,BuyingType,BuyingCostPara,SellingType,SellingCostPara);
%     
% % %     
%     V3cal = reshape((flipud(V3))',NumQ*NumX,1);
%     
%     Ih3 = flipud(reshape(-Aholding*V3cal,NumQ,NumX)');
% 
%     Ib3 = flipud(reshape(-(Abuying*V3cal - Cbuying),NumQ,NumX)');
% 
%     Is3 = flipud(reshape(-(Aselling*V3cal - Cselling),NumQ,NumX)');

% 
% figure 
% hold on
% plot(V1(:,1),'-')
% plot(V2(:,1),'--')
% plot(V3(:,1),'.')
% title({'Solid line is from Markov Chain.';'Dashed line is from solving PDE with Discretization';
%     'Doted line is from solving PDE without Discretization';['NumX = NumQ = ', num2str(NumX)];'q = 0'})
    


%% Holding profit Vs Buying Profit
% figure 
% hold on
% 
% 
% plot(QVec,Ih1(1+13*(2^(k-1)),:))
% plot(QVec,Ib1(1+13*n(k),:),'--')
% 
% title({'Solid line is holding profit and dashed line is buying profit at the price of 0.1';['NumX = ', num2str(NumX); 'NumQ = ', num2str(NumQ)]});
% 
% figure
% hold on 
% 
% plot(QVec(1:(NumQ-1)*0.75+1),Ih1(1+13*i,1:(NumQ-1)*0.75+1))
% plot(QVec(1:(NumQ-1)*0.75+1),Ib1(1+13*i,1:(NumQ-1)*0.75+1),'--')

%% Holding profit minus Buying Profit
% type = ['r','b','k','y'];
% hold on 
% plot(QVec(1:(NumQ-1)*0.75+1),HoldingMinusBuying(1+13*2^(k-1),1:(NumQ-1)*0.75+1),type(k))
% title({'Red line is for NumX = NumQ = 21';'Blue line is for NumX = NumQ = 41';'Black line is for NumX = NumQ = 81';'Log price is 0.1'})%'Yellow line is for NumX = NumQ = 161';
% 


%% Policy plot 
Policy = IMarkov;

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
title({['Optimal policy from policy iteration with beta =', num2str(beta),' difference level =', num2str(DiffLevel)];['NumX = ', num2str(NumX); 'NumQ = ', num2str(NumQ)]});

% plot(QVec,log(CostRatio.*Qmax./QVec),'-')
hold off
% hold on 
% plot(QVec,V1(1+13*(2^(k-1)),:))

end
% V1 = flipud(reshape(V1,NumQ,NumX)');
%% Plot of value function from value iteration with different NumX and NumQ
% Xvec = linspace(Xmin,Xmax,NumX);
% 
% NumP = 1; % The column that we plot
% 
% plot(Xvec,flipud(V0(:,NumP)),'-');
% 
% hold on
% 
% plot(Xvec,flipud(V1(:,NumP)),'--');
% for i = [3,8];
%     NumX = 21 + (i-1)*10;
%     NumQ = NumX;
%     Xvec = linspace(Xmin,Xmax,NumX);
%     [V1,P1] = MarkovChain(alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,NumX,NumQ,MaxIteration,ErrorTol);%V2,P2,
%     V1 = flipud(reshape(V1,NumX,NumQ)');
%     plot(Xvec,flipud(V1(:,NumP)),'--');
% end
% 
% hold off


%% 3d Plot the valuefunction from solvingPDE and value iteration
% Q = Qmin + (Qmax-Qmin)/(NumQ-1) * (repmat((1:NumQ)',1,NumX)-1);
% Q = flipud(Q');
% X = Xmin + (Xmax-Xmin)/(NumX-1) * (repmat((1:NumX),NumQ,1)-1);
% X = flipud(X');
% 
% figure
% mesh(Q,X,V1);
% title('Value function from solving PDE');
% axis([Qmin,Qmax,Xmin,Xmax,0,2000])

% figure %The plot of the most volatile part
% 
% mesh(Q(13:15,:),X(13:15,:),V0(13:15,:))
% figure
% mesh(Q,X,V1);
% title(['Value function from value iteration with Xmin=', num2str(Xmin),' Xmax =', num2str(Xmax)]);
% axis([Qmin,Qmax,Xmin,Xmax,0,2000])

%% Partly fixed Policy iteration 
% [Vpf,Ppf,Npf,Ihpf,Ibpf,Ispf] = MarkovChainPolicyPartlyFixed(alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,NumX,NumQ,MaxIteration,ErrorTol,DiffLevel,CostRatio,I,2);
%  Ipf = reshape(Ppf,NumX,NumQ);
% 
%  Idiff = I - Ipf;
%  
%  % If they are different, display it.
%  if max(max(abs(Idiff))) ~= 0
%      disp(Idiff);
%  end
 
%% Plot of the policy

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
% title({['Optimal policy from policy iteration with beta =', num2str(beta),' difference level =', num2str(DiffLevel)];['NumX = ', num2str(NumX); 'NumQ = ', num2str(NumQ)]});
% 
% % plot(QVec,log(CostRatio.*Qmax./QVec),'-')
% hold off

%% Test the holding value of the value function is symestric or not with respect to alpha

% for k =1 :3
%     figure
%     hold on 
%     plot(linspace(Qmin,Qmax,NumQ),V1(:,19-k),'-');
%     plot(linspace(Qmin,Qmax,NumQ),V1(:,19+k),'--');
%     title({['Solid line is below alpha ', num2str(k*(Xmax-Xmin)/(NumX-1))];['Dashed line is above alpha ', num2str(k*(Xmax-Xmin)/(NumX-1))]})
%     hold off
% end

%% The time that OU process get bigger than the initial time
% Niteration = 1000;
% 
% tau = 0.1;
% T = 100;
% passingtime = 100*ones(Niteration,NumX);
% 
% NumX = 21;
% NumQ = 21;
% 
% XVec = linspace(Xmin,Xmax,NumX);
% 
% 
% for j = 1:NumX
%     x = XVec(j);
%     for i = 1:Niteration
%         X = logPriceGenerator(x,alpha,kappa,sigma,tau,T/tau);
%         t = find(X>x,1,'first')*tau;
%         if isempty(t) == 0
%             passingtime(i,j) = t;
%         end
%     end
% end
%     
% AveragePassingTime = sum(passingtime)/Niteration;
% 
% hold on 
% plot(XVec,AveragePassingTime);
% plot([alpha,alpha],[0,T],'--');
% plotyy(XVec,exp(-beta.*AveragePassingTime),'-.');
% hold off

