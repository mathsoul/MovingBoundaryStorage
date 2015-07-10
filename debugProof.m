InitPara()
global A_hold b_hold A_buy b_buy A_sell b_sell 
[A_hold,b_hold] = GenerateHoldEquation();
[A_buy,b_buy] = GenerateBuyEquation();
[A_sell,b_sell] = GenerateSellEquation();


policy1 = InitMBMPolicy();
policy1(1:(20),1:5,:) = 1;

[value_function1,profit_hold1,profit_buy1,profit_sell1] = SolvePDE3D(policy1);

policy2 = InitMBMPolicy();
policy2(1:(20),1:5,:) = 1;
% policy2(2:end,10:end,:) = 2;

idx = (profit_sell1>-0.0001);
policy2(idx) = 2;

policy = policy2(:,:,1);

[value_function2,~,~,~] = SolvePDE3D(policy2);

diff = value_function1(:,:,1) - value_function2(:,:,1);
