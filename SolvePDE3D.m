function [value_function,profit_hold,profit_buy,profit_sell] = SolvePDE3D(policy)
    
    global NumQ NumX NumS A_hold b_hold A_buy b_buy A_sell b_sell
    
    policy = reshape(policy,[NumQ*NumX*NumS,1]);
        
    [A, b] = GenerateWholeEquation(policy);
        
    value_function_vec = A\b;
    
    value_function = reshape(value_function_vec,NumQ,NumX,NumS);
    
    profit_hold = reshape(A_hold*value_function_vec-b_hold,NumQ,NumX,NumS);

    profit_buy = reshape(-(A_buy*value_function_vec-b_buy),NumQ,NumX,NumS);

    profit_sell = reshape(-(A_sell*value_function_vec-b_sell),NumQ,NumX,NumS);

    profit_buy(NumQ,:,:) = -inf; %When it is full, we can't buy.
    profit_sell(1,:,:) = -inf; %When it is empty, we can't sell.
end












