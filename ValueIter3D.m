% Function ValueIterwSeasonality() computes the optimal policy from value iteration method 
%
% Inputs:
%           Initial value function is set to be all zeros which can be
%           achieved by holding all the time.
%
% Outputs:
%           value_function_VI: value function from value iteration
%           optimal_policy_VI: optimal policy from value iteration
%           n_iter_VI: number of iterations of value iteration method
%           profit_hold_VI:
%           profit_buy_VI:
%           profit_sell:



function [value_function_VI,optimal_policy_VI,n_iter_VI] = ...
    ValueIter3D(start_policy)
    global NumX NumQ NumS MaxIteration ErrorTol
    
    [operator_hold,constant_hold] = GenerateMCHoldOperator();
    [operator_buy,constant_buy] = GenerateMCBuyOperator();
    [operator_sell,constant_sell] = GenerateMCSellOperator();
    
    policy{1} = start_policy;
    
    [operator, constant] = GenerateWholeOperator(policy{1}, ...
    operator_hold, operator_buy, operator_sell, constant_hold, ...
    constant_buy, constant_sell);
    
    policy{1} = reshape(policy{1},[NumQ*NumX*NumS,1]);
    MC{1} = linsolve(eye(NumQ*NumX*NumS) - operator, constant);
   
    converge = 0; %converge or not
    n_iter_VI = 1; %the number of iteration

    while(~converge && n_iter_VI < MaxIteration)

        [MC{2},policy{2}] = max([operator_hold * MC{1} + constant_hold, operator_buy * MC{1} + constant_buy, operator_sell * MC{1} + constant_sell],[],2);
        
        disp(norm(policy{2}-policy{1}))
        
        if(norm(policy{2}-policy{1}) < ErrorTol)
            converge = 1;
        end

        MC{1} = MC{2};
        policy{1} = policy{2};
        
        n_iter_VI = n_iter_VI + 1;
    end

    value_function_VI = MC{1};
    optimal_policy_VI = policy{1} - 1;
    
    value_function_VI = reshape(value_function_VI,[NumQ,NumX,NumS]);
    optimal_policy_VI = reshape(optimal_policy_VI,[NumQ,NumX,NumS]);
end












