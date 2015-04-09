% lambda = buyingCost(Type,QVec,XVec,Qmin,Qmax,CostPara) computes the
% transaction cost at given volume
%Type is a string which indicate the type of the transatction cost. There
%are 3 types which are
% 1. constant; lambda(q) = c_1 + c_2 * e^x; 2. linear; lambda(q) = c_1*(q -
% Qmin + c_2); 3. reciprocal; lambda(q) = c_1 * (Qmax - Qmin + c_2)/( Qmax
% - q + c_2) + c_3 ; Where c_1 c_2 c_3 are constants and given in the
% variable CostPara.



function lambda = buyingCost(Type,QVec,XVec,Qmin,Qmax,CostPara)
    if strcmp(Type,'constant')
        if length(QVec) ~=length(XVec)
            error(1,'QVec and XVec are not of the same length.');
        end
        lambda = CostPara(1) + CostPara(2) .* exp(XVec);
    elseif strcmp(Type,'linear')
        lambda = CostPara(1).* (QVec - Qmin) + CostPara(2);
    elseif strcmp(Type,'reciprocal')
        lambda = CostPara(1) .* (Qmax- Qmin + CostPara(2))./(Qmax - QVec + CostPara(2)) + CostPara(3);
    else
        error(2,'The transaction type is not included.')
    end
    
    if QVec == Qmax
        lambda = inf;
    end
end


 