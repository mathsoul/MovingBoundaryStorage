% q = invBuyingCost calculates the volume when the cost of buying is given.

function q = invSellingCost(cost,Qmax)
    q = Qmax/(cost + 1);
end