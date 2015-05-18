%% This function moves the boundary along X at the same time instead of doing it one by one.
%% MoveIndicator is odd means selling boundary should be moved, otherwise it is even.
function  new_policy = MoveBoundaryAlongX3D(old_policy,profit_buy,profit_sell,which_bd)

    global NumQ NumX NumS
    
    new_policy = old_policy; %initialization
    
    if strcmp(which_bd,'sell')
        for i = 1:NumQ
            for k = 1:NumS
                policy_vec = old_policy(i,:,k);
                profit_sell_vec = profit_sell(i,:,k);
                
                [sell_old_bd,~,~] = getBDFromVec(policy_vec);

                if max(profit_sell_vec) > 0
                    sell_new_bd = find(profit_sell_vec == max(profit_sell_vec),1,'first');
                    
                    new_policy(i,sell_new_bd:(sell_old_bd-1),k) = 2;
                end
            end
        end
    end
    
    if strcmp(which_bd,'buy')
        for i = 1:NumQ
            for k = 1:NumS
                policy_vec = old_policy(i,:,k);
                profit_buy_vec = profit_buy(i,:,k);
                
                [~,buy_old_upper_bd,buy_old_lower_bd] = getBDFromVec(policy_vec);
                if isempty(buy_old_upper_bd)
                    if max(profit_buy_vec) > 0
                        buy_new_upper_bd = ...
                            find(profit_buy_vec == max(profit_buy_vec),1,'last');
                        buy_new_lower_bd = buy_new_upper_bd;
                        new_policy(i,buy_new_upper_bd,k) = 1;
                    end
                else
                    profit_buy_upper_vec = profit_buy(i,buy_old_upper_bd:NumX,k);
                    profit_buy_lower_vec = profit_buy(i,1:buy_old_lower_bd,k);


                    if max(profit_buy_upper_vec) > 0
                        buy_new_upper_bd = buy_old_upper_bd - 1 + ...
                            find(profit_buy_upper_vec == max(profit_buy_upper_vec),1,'last');

                        new_policy(i,(buy_old_upper_bd+1):buy_new_upper_bd,k) = 1;
                    end

                    if max(profit_buy_lower_vec) > 0
                        buy_new_lower_bd = ...
                            find(profit_buy_lower_vec == max(profit_buy_lower_vec),1,'last');

                        new_policy(i,buy_new_lower_bd:(buy_old_upper_bd-1),k) = 1;
                    end
                end
                
                
            end
        end
    end   
end
       
    