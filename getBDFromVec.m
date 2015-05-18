function [sell_bd,buy_upper_bd,buy_lower_bd] = getBDFromVec(policy_vec)
    sell_bd = find(policy_vec == 2,1,'first');
    buy_upper_bd = find(policy_vec == 1,1,'last');
    buy_lower_bd = find(policy_vec == 1,1, 'first');
end