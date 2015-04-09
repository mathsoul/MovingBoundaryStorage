% q = invBuyingCost calculates the volume when the cost of buying is given.

function q = invBuyingCost(cost,Qmax)
    q = Qmax*cost/(cost+1);
end