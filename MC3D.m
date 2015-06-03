InitPara()

% [~,policy_VI,n_iter_VI] = ValueIter3D();

initial_policy = InitMBMPolicy();
% initial_policy = zeros(21,41,12);
[~,policy_PI,n_iter_PI] = PolicyIter3D(initial_policy);

global A_hold b_hold A_buy b_buy A_sell b_sell
[A_hold,b_hold] = GenerateHoldEquation();
[A_buy,b_buy] = GenerateBuyEquation();
[A_sell,b_sell] = GenerateSellEquation();

% value_function_VI = SolvePDE3D(policy_VI);
% value_function_PI = SolvePDE3D(policy_PI);

% save('VI.mat','policy_VI','value_function_VI');
% save('PI.mat','policy_PI','value_function_PI');
n_PI = length(policy_PI);

for i = 1:n_PI
    value_function_PI{i} = SolvePDE3D(policy_PI{i});
end
save('PI_policy_seq','policy_PI','value_function_PI');
