InitPara()

[~,policy_VI,n_iter_VI] = ValueIter3D();

[~,policy_PI,n_iter_PI] = PolicyIter3D(policy_VI);

global A_hold b_hold A_buy b_buy A_sell b_sell
[A_hold,b_hold] = GenerateHoldEquation();
[A_buy,b_buy] = GenerateBuyEquation();
[A_sell,b_sell] = GenerateSellEquation();

value_function_VI = SolvePDE3D(policy_VI);
value_function_PI = SolvePDE3D(policy_PI);

save('VI.mat','policy_VI','value_function_VI');
save('PI.mat','policy_PI','value_function_PI');
