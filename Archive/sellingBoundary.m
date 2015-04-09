% This one is wrong because the assumption that the boundary can be
% analytically calculated is wrong.
% x = sellingBoundary() computes the selling boundary volume q given log price x.

function q = sellingBoundary(x,alpha,beta,kappa,sigma,Qmax)
    
    %original boundary

%      cost = -(0.5*sigma^2*exp(x) + kappa*(alpha-x)*exp(x) - beta*exp(x))/(beta);


    %new boundary
    x = x-0.2;
    cost = -(0.5*sigma^2*exp(x) + kappa*(alpha-x)*exp(x) - beta*exp(x))/(4*beta);
    if cost < 0
        q = 100;
        return 
    end
    q = invSellingCost(cost,Qmax);
end
