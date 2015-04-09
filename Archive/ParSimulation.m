% [Valuefunction,StandardDeviation,BuyingTimes,SellingTimes,NegativeSellingTimes,TotalBuyingCosts,TotalSellingCosts]
% = Simulation(x,q,alpha,beta,kappa,sigma,Qmax,Qmin) Computes the value
% function via doing simulation. Here the boundary is fixed and analytical
% which is discribed in another two functions buyingBoundary and
% sellingBoundary.
%
% The inputs are:
%   x is the initial log price
%   q is the initial volume
%   alpha is the long term mean value
%   beta is the discount factor
%   kappa is the rate the value reverts to the mean
%   sigma is the volatility 
%   Qmax is the maximum value of volume
%   Qmin is the minimum value of volume
%   T is the time limit of the simulatio
%   N is the number of simulation 
%   tau is the increasement of time
%
% and returns:
%   Valuefunction is the value of vaule function begining with log price x
%   and volume q. 
%   StandardDeviation is the standard deviation of the N times simulation.
%   BuyingTimes is the average times of buying.
%   SellingTimes is the average times of selling.
%   NegativeSellingTimes is the average time of selling with negative
%   profit.
%   TotalBuyingCosts is the average value of buying costs.
%   TotalSellingCosts is the average value of selling costs.
function [Valuefunction,StandardDeviation,BuyingTimes,SellingTimes,...
    NegativeSellingTimes,TotalBuyingCosts,TotalSellingCosts] = ...
    ParSimulation(x,q,alpha,beta,kappa,sigma,Qmax,Qmin,T,N,tau)
%% SIMULATION

% matlabpool(4)

% Initialization
V = zeros(N,1); %value function matrix

[q,Valuefunction,BuyingTimes,SellingTimes,...
    NegativeSellingTimes,TotalBuyingCosts,TotalSellingCosts] ...
    =InitialTransaction(x,q,alpha,beta,kappa,sigma,Qmax,Qmin,T,N,tau);

V(:) = Valuefunction;

parfor n = 1:N
    X = logPriceGenerator(x,alpha,kappa,sigma,tau,T/tau); %log price vector w.r.t time
    Q = zeros(T/tau,1); %volume vector w.r.t. time
    Q(1) = q;
    for k = 1:(T/tau-1)
%         Buy, sell or hold
        if( Q(k) <  buyingBoundary(X(k+1),alpha,beta,kappa,sigma,Qmax)) % buy
            BuyingTimes = BuyingTimes +1;
            Q(k+1) =buyingBoundary(X(k+1),alpha,beta,kappa,sigma,Qmax);
            cost = Qmax* log((Qmax-Q(k))/(Qmax-Q(k+1))) - (Q(k+1) - Q(k));
            TotalBuyingCosts = TotalBuyingCosts + cost *exp(-beta*(k+1)*tau) ;
            V(n) =  V(n) - (exp(X(k+1))*(Q(k+1)-Q(k)) + cost)*exp(-beta*(k+1)*tau);

        elseif(Q(k)> sellingBoundary(X(k+1),alpha,beta,kappa,sigma,Qmax)) %Sell
            SellingTimes = SellingTimes +1;
            Q(k+1) = sellingBoundary(X(k+1),alpha,beta,kappa,sigma,Qmax);
            cost = Qmax * log(Q(k)/Q(k+1))-(Q(k)-Q(k+1));
            TotalSellingCosts = TotalSellingCosts + cost*exp(-beta*(k+1)*tau);
            V(n) = V(n)- (exp(X(k+1))*(Q(k+1)-Q(k)) + cost)*exp(-beta*(k+1)*tau);
            if(exp(X(k+1))<sellingCost(Q(k+1),Qmax))
                NegativeSellingTimes = NegativeSellingTimes+1;
            end
        else
            V(n) = V(n);
            Q(k+1) = Q(k);
        end
    end
%     plot(tau:tau:T,X)
%     hold on
%     scatter(Q,X)
    
end

BuyingTimes = BuyingTimes/N;

SellingTimes = SellingTimes/N;

NegativeSellingTimes = NegativeSellingTimes/N;

TotalBuyingCosts = TotalBuyingCosts/N;

TotalSellingCosts = TotalSellingCosts/N;

StandardDeviation = sqrt(var(V)/N);

Valuefunction = mean(V);

% matlabpool close

end











