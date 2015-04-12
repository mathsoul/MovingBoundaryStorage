% Function ValueIterwSeasonality() computes the optimal policy from value iteration method 
%
% Inputs:
%           initial_region is set to be all holding
%
% Outputs:
%           value_function_VI: value function from value iteration
%           optimal_policy_VI: optimal policy from value iteration
%           n_iter_VI: number of iterations of value iteration method
%           profit_hold_VI:
%           profit_buy_VI:
%           profit_sell:



function [value_function_VI,optimal_policy_VI,n_iter_VI] = ...
    ValueIterwSeasonality()
    global NumX NumQ NumS MaxIteration ErrorTol
    
    [operator_hold,constant_hold] = GenerateMCHoldOperator();
    [operator_buy,constant_buy] = GenerateMCBuyOperator();
    [operator_sell,constant_sell] = GenerateMCSellOperator();


    MC{1}= zeros(NumQ*NumX*NumS,1);
    converge = 1; %converge or not
    n_iter_VI = 1; %the number of iteration

    while(~converge && n_iter_VI<MaxIteration)

        [MC{2},Policy] = max([operator_hold * MC{1} + constant_hold, operator_buy * MC{1} + constant_buy, operator_sell * MC{1} + constant_sell],[],2);

        if(norm(MC{2}-MC{1}) < ErrorTol)
            converge = 1;
        end

        MC{1} = MC{2};

        n_iter_VI = n_iter_VI + 1;
    end

    value_function_VI = MC{1};
    optimal_policy_VI = Policy - 1;

end












