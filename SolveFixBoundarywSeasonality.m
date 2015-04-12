% Function SolvePDEwSeasonality() solve the fixed boundary problem with
% seasonality
% Inputs: 
%           policy: 2 sell, 1 buy, 0 hold
% Outputs:
%           value_function
%           profit_hold
%           profit_sell
%           profit_buy


function [value_function,profit_hold,profit_sell,profit_buy] = SolveFixBoundarywSeasonality(policy)
    
    [A_hold,b_hold] = GenerateHoldEquation();
    [A_buy,b_buy] = GenerateBuyEquation();
    [A_sell,b_sell] = GenerateSellEquation();

    
    [A,b] = GenerateWholeEquation(policy, A_hold, A_buy, A_sell, b_hold, b_buy, b_sell);
    
    value_functionVec = linsolve(A,b);
    value_function = Reshape4Disp(value_functionVec);
    profit_hold = Reshape4Disp(-(A_hold*value_functionVec - b_hold));
    profit_buy = Reshape4Disp(-(A_buy*value_functionVec - b_buy));
    profit_sell = Reshape4Disp(-(A_sell*value_functionVec - b_sell));
end




    



               