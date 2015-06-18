InitPara()

global A_hold b_hold A_buy b_buy A_sell b_sell AbsDiff
[A_hold,b_hold] = GenerateHoldEquation();
[A_buy,b_buy] = GenerateBuyEquation();
[A_sell,b_sell] = GenerateSellEquation();

%% Move the boundary along X at the same time
initial_policy = InitMBMPolicy();

old_policy = initial_policy;

[value_function{1},profit_hold{1},profit_buy{1},profit_sell{1}] = SolvePDE3D(old_policy);  

clear new_policy

new_policy{1} = MoveBoundaryAlongX3D(old_policy,profit_buy{1},profit_sell{1},'sim');

counter = 1;

while norm3D(new_policy{counter}-old_policy)>0
   old_policy = new_policy{counter};
   counter = counter +1;
   [value_function{counter},profit_hold{counter},profit_buy{counter},...
        profit_sell{counter}] = SolvePDE3D(old_policy);
   new_policy{counter} = ...
       MoveBoundaryAlongX3D(old_policy,profit_buy{counter},profit_sell{counter},'sim');
end

clear policy_MBM_sim

policy_MBM_sim = new_policy{counter};   
value_function_MBM_sim = value_function{counter};
save(sprintf('MBM_sim_%.1f.mat', AbsDiff),'policy_MBM_sim','value_function_MBM_sim');

policy_MBM_sim = new_policy;
value_function_MBM_sim = value_function;
save(sprintf('MBM_sim_policy_seq%.1f.mat', AbsDiff),'policy_MBM_sim','value_function_MBM_sim');

% for j = 1:length(new_policy)
%     figure
%     BoundaryLimit(new_policy{15}(:,:,11))
%     hold on
%     plot([35,35],[-100,100])
%     hold off
% end

