% mu = sellingCost(Type,QVec,XVec,Qmin,Qmax,CostPara)
%computes the selling cost at given volume
%Type is a string which indicate the type of the transatction cost.
%There are 3 types which are
% 1. constant; mu(q) = c_1 + c_2 * e^x;
% 2. linear; mu(q) = c_1*(Qmax - q + c_2);
% 3. reciprocal; mu(q) = c_1 * (Qmax - Qmin + c_2)/( q - Qmin + c_2) + c_3 ;
% Where c_1 c_2 c_3 are constants and given in the variable CostPara.
function mu = sellingCost(Type,QVec,XVec,Qmin,Qmax,CostPara)
    if strcmp(Type,'constant')
        if length(QVec) ~=length(XVec)
            error(1,'QVec and XVec are not of the same length.');
        end
        mu = CostPara(1) + CostPara(2) .* exp(XVec);
    elseif strcmp(Type,'linear')
        mu = CostPara(1).* (Qmax - QVec)  + CostPara(2);
    elseif strcmp(Type,'reciprocal')
        mu = CostPara(1) .* (Qmax-Qmin+CostPara(2))./(QVec - Qmin + CostPara(2)) + CostPara(3);
    else
        error(2,'The transaction type is not included.')
    end
    
    if QVec == Qmin
        mu = inf;
    end
end