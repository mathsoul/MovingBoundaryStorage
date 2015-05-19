%cost = AggregatedSellingCost(Type,qlower,qupper,x,Qmax,Qmin,CostPara)
%Type is a string which indicate the type of the transatction cost.
%There are 3 types which are
% 1. constant; mu(q) = c_1 + c_2 * e^x;
% 2. linear; mu(q) = c_1*(Qmax - q) + c_2 ;
% 3. reciprocal; mu(q) = c_1 * (Qmax - Qmin + c_2)/( q - Qmin + c_2) + c_3 ;
% Where c_1 c_2 c_3 are constants and given in the variable CostPara.


function cost = AggregatedSellingCost(Type,qlower,qupper,x,Qmax,Qmin,NumQ,CostPara)
    if strcmp(Type,'constant')
        cost = (CostPara(1) + CostPara(2) * exp(x))*(qupper - qlower);
    elseif strcmp(Type,'linear')
        cost = CostPara(1)* ( -(qupper^2 - qlower^2)/2 + Qmax *(qupper - qlower)) + CostPara(2)*(qupper - qlower);
    elseif strcmp(Type,'reciprocal')
        if CostPara(2) == 0 && qlower == Qmin % To avoid infinite transaction cost
            dq = (Qmax - Qmin)/(NumQ - 1);
            cost = AggregatedSellingCost(Type,Qmin+dq,qupper,x,Qmax,Qmin,NumQ,CostPara) + dq * (CostPara(1) * (Qmax - Qmin)/dq + CostPara(3));
            cost = cost*exp(x);
        else
            cost = CostPara(1) * (Qmax-Qmin+CostPara(2)) * log((qupper - Qmin +CostPara(2)) /(qlower - Qmin + CostPara(2))) + CostPara(3)*(qupper - qlower);
            cost = cost*exp(x);
        end
    else
        error(2,'The transaction type is not included.')
    end
end