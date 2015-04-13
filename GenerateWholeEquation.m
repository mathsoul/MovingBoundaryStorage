% Function () generates the linear equations based on the region's behavior
% Inputs:
%           A_hold (A_sell, A_buy) is the A matrix of holding (selling,
%           buying)
%           b_hold (b_sell, b_buy) is the b vec of holding (selling,
%           buying)
% Outputs:
%           A is the A matrix of the whole system
%           b is the b vec of the whole system


function [A,b] = GenerateWholeEquation(policy, A_hold, A_buy, A_sell, b_hold, b_buy, b_sell)
    
    global NumQ NumX NumS
    
    [A, b] = InitEquation();
    
    for ijk=1:NumQ*NumX*NumS

        if(policy(ijk) == 1)
           A(ijk,:) = A_buy(ijk,:);
           b(ijk) = b_buy(ijk);
        elseif(policy(ijk) == 2)
            A(ijk,:) = A_sell(ijk,:);
            b(ijk) = b_sell(ijk);
        else
            A(ijk,:) = A_hold(ijk,:);
            b(ijk) = b_hold(ijk);
        end
    end
end