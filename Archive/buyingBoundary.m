% This one is wrong because the assumption that the boundary can be
% analytically calculated is wrong.
% x = buyingBoundary() computes the buying boundary volume q given log price x.

function q = buyingBoundary(x,alpha,beta,kappa,sigma,Qmax)
    %original boundary
%      cost = (0.5*sigma^2*exp(x) + kappa*(alpha-x)*exp(x) - beta*exp(x))/(beta);


    %new boundary
    x = x+0.2;
    cost = (0.5*sigma^2*exp(x) + kappa*(alpha-x)*exp(x) - beta*exp(x))/(10*beta);
    if cost < 0
        q = 0;
        return
    end
    q = invBuyingCost(cost,Qmax);
end
