%

function Q = VolumeUnderGivenPolicy(x,q,alpha,beta,kappa,sigma,Qmax,Qmin,T,tau)
    X = logPriceGenerator(x,alpha,kappa,sigma,T,tau);
    Q = zeros(T/tau,1); %volume vector w.r.t. time
    Q(1) = q;
        if( Q(1) <  buyingBoundary(X(1),alpha,beta,kappa,sigma,Qmax)) % buy at initial time
            BuyingTimes = BuyingTimes+1;
            Q(1) =buyingBoundary(X(1),alpha,beta,kappa,sigma,Qmax);
            cost = Qmax* log((Qmax-q)/(Qmax-Q(1))) - (Q(1) - q);
            TotalBuyingCosts = TotalBuyingCosts + cost  ;
            V(n,1) =  - (exp(X(1))*(Q(1)-q) + cost);

        elseif(Q(1)> sellingBoundary(X(1),alpha,beta,kappa,sigma,Qmax)) %Sell at initial time
            SellingTimes = SellingTimes +1;
            Q(1) = sellingBoundary(X(1),alpha,beta,kappa,sigma,Qmax);
            cost = Qmax * log(q/Q(1))-(q-Q(1));
            TotalSellingCosts = TotalSellingCosts + cost;
            V(n,1) =  -(exp(X(1))*(Q(1)-q) + cost);
            if(exp(X(1))<sellingCost(Q(1),Qmax))
                NegativeSellingTimes = NegativeSellingTimes+1;
            end
        end
    for k = 1:(T/tau-1)
%         Buy, sell or hold
        if( Q(k) <  buyingBoundary(X(k+1),alpha,beta,kappa,sigma,Qmax)) % buy
            BuyingTimes = BuyingTimes +1;
            Q(k+1) =buyingBoundary(X(k+1),alpha,beta,kappa,sigma,Qmax);
            cost = Qmax* log((Qmax-Q(k))/(Qmax-Q(k+1))) - (Q(k+1) - Q(k));
            TotalBuyingCosts = TotalBuyingCosts + cost *exp(-beta*(k+1)*tau) ;
            V(n,k + 1) =  V(n,k) - (exp(X(k+1))*(Q(k+1)-Q(k)) + cost)*exp(-beta*(k+1)*tau);

        elseif(Q(k)> sellingBoundary(X(k+1),alpha,beta,kappa,sigma,Qmax)) %Sell
            SellingTimes = SellingTimes +1;
            Q(k+1) = sellingBoundary(X(k+1),alpha,beta,kappa,sigma,Qmax);
            cost = Qmax * log(Q(k)/Q(k+1))-(Q(k)-Q(k+1));
            TotalSellingCosts = TotalSellingCosts + cost*exp(-beta*(k+1)*tau);
            V(n,k+1) = V(n,k)- (exp(X(k+1))*(Q(k+1)-Q(k)) + cost)*exp(-beta*(k+1)*tau);
            if(exp(X(k+1))<sellingCost(Q(k+1),Qmax))
                NegativeSellingTimes = NegativeSellingTimes+1;
            end
        else
            V(n,k+1) = V(n,k);
            Q(k+1) = Q(k);
        end
    end
end