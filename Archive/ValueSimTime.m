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
function Valuefunction = ...
    ValueSimTime(x,q,alpha,beta,kappa,sigma,Qmax,Qmin,T,N,tau)
%% SIMULATION

% Initialization
V = zeros(N,T/tau); %value function matrix

[q,InitialValue,~]=InitialTransaction(x,q,alpha,beta,kappa,sigma,Qmax,Qmin,T,N,tau);

V(:,1) = InitialValue;

for n = 1:N
    X = logPriceGenerator(x,alpha,kappa,sigma,tau,T/tau); %log price vector w.r.t time
    X(1) = x;
    Q = zeros(T/tau,1); %volume vector w.r.t. time
    Q(1) = q;
    for k = 1:(T/tau-1)
%         Buy, sell or hold
        if( Q(k) <  buyingBoundary(X(k+1),alpha,beta,kappa,sigma,Qmax)) % buy
            Q(k+1) =buyingBoundary(X(k+1),alpha,beta,kappa,sigma,Qmax);
            cost = Qmax* log((Qmax-Q(k))/(Qmax-Q(k+1))) - (Q(k+1) - Q(k));
            V(n,k + 1) =  V(n,k) - (exp(X(k+1))*(Q(k+1)-Q(k)) + cost)*exp(-beta*(k+1)*tau);

        elseif(Q(k)> sellingBoundary(X(k+1),alpha,beta,kappa,sigma,Qmax)) %Sell
            Q(k+1) = sellingBoundary(X(k+1),alpha,beta,kappa,sigma,Qmax);
            cost = Qmax * log(Q(k)/Q(k+1))-(Q(k)-Q(k+1));
            V(n,k+1) = V(n,k)- (exp(X(k+1))*(Q(k+1)-Q(k)) + cost)*exp(-beta*(k+1)*tau);
        else
            V(n,k+1) = V(n,k);
            Q(k+1) = Q(k);
        end
    end
end

Valuefunction = mean(V);

end











