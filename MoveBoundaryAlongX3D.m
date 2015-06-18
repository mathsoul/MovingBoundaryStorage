%% This function moves the boundary along X at the same time instead of doing it one by one.
%% MoveIndicator is odd means selling boundary should be moved, otherwise it is even.
function  new_policy = MoveBoundaryAlongX3D(old_policy,profit_buy,profit_sell,which_bd)

    global NumQ NumX NumS SellAbsDiff BuyAbsDiff
    
    new_policy = old_policy; %initialization
    
    if strcmp(which_bd,'sell')
        for i = 1:NumQ
            for k = 1:NumS
                policy_vec = old_policy(i,:,k);
                [sell_old_bd,~,~] = getBDFromVec(policy_vec);
                
                profit_sell_vec = profit_sell(i,1:(sell_old_bd-1),k);
                % select sequentially positive sub vector
                last_neg_index = 0;
                if min(profit_sell_vec)<0
                    last_neg_index = find(profit_sell_vec < 0,1,'last');
                    profit_sell_vec = profit_sell_vec((last_neg_index+1):length(profit_sell_vec));
                end
                
                if max(profit_sell_vec) > SellAbsDiff
                    sell_new_bd = find(profit_sell_vec == max(profit_sell_vec),1,'first') + last_neg_index;
                    
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
                    if max(profit_buy_vec) > BuyAbsDiff
                        buy_new_upper_bd = ...
                            find(profit_buy_vec == max(profit_buy_vec),1,'last');
                        buy_new_lower_bd = buy_new_upper_bd;
                        new_policy(i,buy_new_upper_bd,k) = 1;
                    end
                else
                    profit_buy_upper_vec = profit_buy(i,(buy_old_upper_bd+1):NumX,k);
                    % select sequentially positive sub vector
                    if min(profit_buy_upper_vec)<0
                        last_neg_index = find(profit_buy_upper_vec<0,1,'last');
                        profit_buy_upper_vec = ...
                            profit_buy_upper_vec(1:(last_neg_index-1));
                    end

                   if max(profit_buy_upper_vec) > BuyAbsDiff
                        buy_new_upper_bd = buy_old_upper_bd  + ...
                            find(profit_buy_upper_vec == max(profit_buy_upper_vec),1,'last');

                        new_policy(i,(buy_old_upper_bd+1):buy_new_upper_bd,k) = 1;
                    end
                    
                    profit_buy_lower_vec = profit_buy(i,1:(buy_old_lower_bd-1),k);
                    % select sequentially positive sub vector
                    if min(profit_buy_lower_vec)<0
                        first_neg_index = find(profit_buy_lower_vec<0,1,'first');
                        profit_buy_lower_vec = ...
                            profit_buy_lower_vec(1:(first_neg_index-1));
                    end
 
                    if max(profit_buy_lower_vec) > BuyAbsDiff
                        buy_new_lower_bd = ...
                            find(profit_buy_lower_vec == max(profit_buy_lower_vec),1,'last');

                        new_policy(i,buy_new_lower_bd:(buy_old_upper_bd-1),k) = 1;
                    end
                end
            end
        end
    end
    
    if strcmp(which_bd,'sim')
        for i = 1:NumQ
            for k = 1:NumS
                % move the selling boundary
                policy_vec = old_policy(i,:,k);
                [sell_old_bd,~,~] = getBDFromVec(policy_vec);
                
                profit_sell_vec = profit_sell(i,1:(sell_old_bd-1),k);
                % select sequentially positive sub vector
                last_neg_index = 0;
                if min(profit_sell_vec)<0
                    last_neg_index = find(profit_sell_vec < 0,1,'last');
                    profit_sell_vec = profit_sell_vec((last_neg_index+1):length(profit_sell_vec));
                end
                
                if max(profit_sell_vec) > SellAbsDiff
                    sell_new_bd = find(profit_sell_vec == max(profit_sell_vec),1,'first') + last_neg_index;
                    
                    new_policy(i,sell_new_bd:(sell_old_bd-1),k) = 2;
                end
                
                % move the buying boundary
                policy_vec = old_policy(i,:,k);
                profit_buy_vec = profit_buy(i,:,k);
                
                [~,buy_old_upper_bd,buy_old_lower_bd] = getBDFromVec(policy_vec);
                if isempty(buy_old_upper_bd)
                    if max(profit_buy_vec) > BuyAbsDiff
                        buy_new_upper_bd = ...
                            find(profit_buy_vec == max(profit_buy_vec),1,'last');
                        buy_new_lower_bd = buy_new_upper_bd;
                        new_policy(i,buy_new_upper_bd,k) = 1;
                    end
                else
                    profit_buy_upper_vec = profit_buy(i,(buy_old_upper_bd+1):NumX,k);
                    % select sequentially positive sub vector
                    if min(profit_buy_upper_vec)<0
                        last_neg_index = find(profit_buy_upper_vec<0,1,'last');
                        profit_buy_upper_vec = ...
                            profit_buy_upper_vec(1:(last_neg_index-1));
                    end

                   if max(profit_buy_upper_vec) > BuyAbsDiff
                        buy_new_upper_bd = buy_old_upper_bd  + ...
                            find(profit_buy_upper_vec == max(profit_buy_upper_vec),1,'last');

                        new_policy(i,(buy_old_upper_bd+1):buy_new_upper_bd,k) = 1;
                    end
                    
                    profit_buy_lower_vec = profit_buy(i,1:(buy_old_lower_bd-1),k);
                    % select sequentially positive sub vector
                    if min(profit_buy_lower_vec)<0
                        first_neg_index = find(profit_buy_lower_vec<0,1,'first');
                        profit_buy_lower_vec = ...
                            profit_buy_lower_vec(1:(first_neg_index-1));
                    end
 
                    if max(profit_buy_lower_vec) > BuyAbsDiff
                        buy_new_lower_bd = ...
                            find(profit_buy_lower_vec == max(profit_buy_lower_vec),1,'last');

                        new_policy(i,buy_new_lower_bd:(buy_old_upper_bd-1),k) = 1;
                    end
                end
            end
        end
    end
end
       
    