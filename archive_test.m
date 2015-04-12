%Policy iteration
Policy{1} = PolicyValueIteration;
converge = 1; %converge or not
n_iter_VI = 1; %the number of iteration

while(converge && n_iter_VI<MaxIteration)
    if(Policy{n_iter_VI}(ijk) == 1)
       A(ijk,:) = A_buy(ijk,:);
       b(ijk) = b_buy(ijk);
    elseif(Policy{n_iter_VI}(ijk) == 2)
        A(ijk,:) = A_sell(ijk,:);
        b(ijk) = b_sell(ijk);
    else
        A(ijk,:)=A_hold(ijk,:);
        b(ijk) = b_hold(ijk);
    end
    MC{n_iter_VI} = linsolve(eye(NumQ*NumX*NumS) - A,b); % I am not sure what this should be
    V_hold = A_hold * MC{n_iter_VI} + b_hold;
    V_buy = A_buy * MC{n_iter_VI} + b_buy;
    V_sell = A_sell * MC{n_iter_VI} + b_sell;
    [a,b] = max([V_hold,V_buy,V_sell],[],2);
    
    
    
    for ijk = 1:NumQ*NumX*NumS
        if(a(ijk)> ( 1+ DiffLevel)*MC{n_iter_VI}(ijk))
            Policy{n_iter_VI+1}(ijk,1) = b(ijk) - 1;
        else
            Policy{n_iter_VI+1}(ijk,1) = Policy{n_iter_VI}(ijk,1);
        end
    end

    if(norm(Policy{n_iter_VI+1}-Policy{n_iter_VI}) < ErrorTol)
        converge = 0;
    end
    n_iter_VI = n_iter_VI+1;
end

value_function_VI = MC{n_iter_VI-1};
optimal_policy_VI = Policy{n_iter_VI-1};
n_iter_VI = n_iter_VI-1;

profit_hold_VI = reshape4disp(V_hold,NumQ,NumX,NumS);

profit_buy_VI = reshape4disp(V_hold,NumQ,NumX,NumS);

profit_sell_VI = reshape4disp(V_hold,NumQ,NumX,NumS);