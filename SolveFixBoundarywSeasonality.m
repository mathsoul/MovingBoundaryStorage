% Function SolvePDEwSeasonality() solve the fixed boundary problem with
% seasonality
% Inputs: 
%           region: the indicator of policies
% Outputs:
%           value_function
%           profit_hold
%           profit_sell
%           profit_buy


function [value_function,profit_hold,profit_sell,profit_buy] = SolveFixBoundarywSeasonality(region)
    
    global NumX NumQ NumS
    
    [~, indexVecQ, indexVecX, indexVecS] = NodeIndex();
    
    [A, b] = InitEquation();
    
    [A_hold,b_hold] = GenerateHoldEquation();
    [A_sell,b_sell] = GenerateSellEquation();
    [A_buy,b_buy] = GenerateBuyEquation();
    
    for ijk=1:NumQ*NumX*NumS
        i = indexVecQ(ijk);
        j = indexVecX(ijk);
        k = indexVecS(ijk);

        if(region(i,j,k) == 1)
           A(ijk,:) = A_buy(ijk,:);
           b(ijk) = b_buy(ijk);
        elseif(region(i,j,k) == 2)
            A(ijk,:) = A_sell(ijk,:);
            b(ijk) = b_sell(ijk);
        else
            A(ijk,:)=A_hold(ijk,:);
            b(ijk) = b_hold(ijk);
        end
    end
    
    value_functionVec = linsolve(A,b);
    value_function = Reshape4Disp(value_functionVec);
    profit_hold = Reshape4Disp(-(A_hold*value_functionVec - b_hold));
    profit_buy = Reshape4Disp(-(A_buy*value_functionVec - b_buy));
    profit_sell = Reshape4Disp(-(A_sell*value_functionVec - b_sell));
end




    



               