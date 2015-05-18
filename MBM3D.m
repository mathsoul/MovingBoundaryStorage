InitPara()

global A_hold b_hold A_buy b_buy A_sell b_sell AbsDiff
[A_hold,b_hold] = GenerateHoldEquation();
[A_buy,b_buy] = GenerateBuyEquation();
[A_sell,b_sell] = GenerateSellEquation();

%% Move the boundary along X at the same time
initial_policy = InitMBMPolicy();

old_policy = initial_policy;

[value_function{1},profit_hold{1},profit_buy{1},profit_sell{1}] = SolvePDE3D(old_policy);  

new_policy{1} = MoveBoundaryAlongX3D(old_policy,profit_buy{1},profit_sell{1},'sell');

counter = 1;

while norm3D(new_policy{counter}-old_policy)>0
   % Move buying boundary
   old_policy = new_policy{counter};
   counter = counter +1;
   [value_function{counter},profit_hold{counter},profit_buy{counter},...
        profit_sell{counter}] = SolvePDE3D(old_policy);
   new_policy{counter} = ...
       MoveBoundaryAlongX3D(old_policy,profit_buy{counter},profit_sell{counter},'buy');
   
   % Move selling boundary
   old_policy = new_policy{counter};
   counter = counter +1;
   [value_function{counter},profit_hold{counter},profit_buy{counter},...
        profit_sell{counter}] = SolvePDE3D(old_policy);
   new_policy{counter} = ...
       MoveBoundaryAlongX3D(old_policy,profit_buy{counter},profit_sell{counter},'sell');
   
   % Recover previous policy
   old_policy = new_policy{counter -2};
end


policy_MBM = new_policy{counter};   
value_function_MBM = value_function{counter};
save(sprintf('MBM_%.1f.mat', AbsDiff),'policy_MBM','value_function_MBM');


% for j = 1:length(new_policy)
%     figure
%     BoundaryLimit(new_policy{15}(:,:,11))
%     hold on
%     plot([35,35],[-100,100])
%     hold off
% end

