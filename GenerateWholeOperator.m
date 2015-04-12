% Function () generates the linear equations based on the policy
% Inputs:
%           A_hold (A_sell, A_buy) is the A matrix of holding (selling,
%           buying)
%           b_hold (b_sell, b_buy) is the b vec of holding (selling,
%           buying)
% Outputs:
%           A is the A matrix of the whole system
%           b is the b vec of the whole system


function [operator,constant] = GenerateWholeOperator(policy, operator_hold, operator_buy, operator_sell, constant_hold, constant_buy, constant_sell)
    
    global NumQ NumX NumS
    [~, indexVecQ, indexVecX, indexVecS] = NodeIndex();
    
    [operator, constant] = InitOperator();
    
    for ijk=1:NumQ*NumX*NumS
        i = indexVecQ(ijk);
        j = indexVecX(ijk);
        k = indexVecS(ijk);

        if(policy(i,j,k) == 1)
           operator(ijk,:) = operator_buy(ijk,:);
           constant(ijk) = constant_buy(ijk);
        elseif(policy(i,j,k) == 2)
            operator(ijk,:) = operator_sell(ijk,:);
            constant(ijk) = constant_sell(ijk);
        else
            operator(ijk,:)=operator_hold(ijk,:);
            constant(ijk) = constant_hold(ijk);
        end
    end
    
end