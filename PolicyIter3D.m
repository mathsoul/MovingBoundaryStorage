% Function PolicyIterwSeasonality() computes the optimal policy from policy
% iteration method with a start policy
%
% Inputs:
%           start_policy
% Outputs:
%           value_function_PI: value function from value iteration
%           optimal_policy_PI: optimal policy from value iteration
%           n_iter_PI: number of iterations of value iteration method
%           profit_hold_VI:
%           profit_buy_VI:
%           profit_sell:



function [value_function_PI,optimal_policy_PI,n_iter_PI] = ...
    PolicyIter3D(start_policy)
    global NumQ NumX NumS MaxIteration ErrorTol DiffLevel
    
    [operator_hold,constant_hold] = GenerateMCHoldOperator();
    [operator_buy,constant_buy] = GenerateMCBuyOperator();
    [operator_sell,constant_sell] = GenerateMCSellOperator();
    
    policy{1} = reshape(start_policy,[NumQ*NumX*NumS,1]);
    converge = 0; %converge or not
    n_iter_PI = 1; %the number of iteration

    while(~converge && n_iter_PI < MaxIteration)
        
        [operator, constant] = GenerateWholeOperator(policy{n_iter_PI}, ...
            operator_hold, operator_buy, operator_sell, constant_hold, ...
            constant_buy, constant_sell);
        
        MC = linsolve(eye(NumQ*NumX*NumS) - operator, constant);
        
        [max_value,which_max] = max([operator_hold * MC + constant_hold, operator_buy * MC + constant_buy, operator_sell * MC + constant_sell],[],2);
        
        for ijk = 1:NumQ*NumX*NumS
            if( max_value(ijk) > ( 1+ DiffLevel) * MC(ijk))
                policy{n_iter_PI+1}(ijk,1) = which_max(ijk) - 1;
            else
                policy{n_iter_PI+1}(ijk,1) = policy{n_iter_PI}(ijk,1);
            end
        end
        
        
        if(norm(policy{n_iter_PI + 1} - policy{n_iter_PI}) < ErrorTol)
            converge = 1;
        end

        n_iter_PI = n_iter_PI + 1;
    end

    value_function_PI = MC;
%     optimal_policy_PI = policy{n_iter_PI};
    optimal_policy_PI = policy;
    for i = 1:n_iter_PI
        optimal_policy_PI{i} = reshape(optimal_policy_PI{i},[NumQ,NumX,NumS]);
    end
    
    value_function_PI = reshape(value_function_PI,[NumQ,NumX,NumS]);
%     optimal_policy_PI = reshape(optimal_policy_PI,[NumQ,NumX,NumS]);
end












