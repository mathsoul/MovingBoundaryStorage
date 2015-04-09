
function [q,Valuefunction,BuyingTimes,SellingTimes,...
    NegativeSellingTimes,InitialBuyingCosts,InitialSellingCosts] = ...
    InitialTransaction(xInitial,qInitial,alpha,beta,kappa,sigma,Qmax,Qmin,T,N,tau)
    Valuefunction = 0; 
    BuyingTimes = 0;
    SellingTimes = 0;
    NegativeSellingTimes = 0;
    InitialBuyingCosts = 0;
    InitialSellingCosts = 0;

    if( qInitial <  buyingBoundary(xInitial,alpha,beta,kappa,sigma,Qmax)) % buy at initial time
        BuyingTimes = 1;
        q = buyingBoundary(xInitial,alpha,beta,kappa,sigma,Qmax);
        cost = Qmax* log((Qmax-qInitial)/(Qmax-q)) - (q - qInitial);
        InitialBuyingCosts = cost  ;
        Valuefunction=  - (exp(xInitial)*(q-qInitial) + cost);
    elseif(qInitial> sellingBoundary(xInitial,alpha,beta,kappa,sigma,Qmax)) %Sell at initial time
        SellingTimes = 1;
        q = sellingBoundary(xInitial,alpha,beta,kappa,sigma,Qmax);
        cost = Qmax * log(qInitial/q)-(qInitial-q);
        InitialSellingCosts =  cost;
        Valuefunction(:) =  -(exp(xInitial)*(q-qInitial) + cost);
        if(exp(xInitial)<sellingCost(q,Qmax))
            NegativeSellingTimes = NegativeSellingTimes+1;
        end
    else
        q = qInitial;
    end
    

end