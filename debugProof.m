policy1 = InitMBMPolicy();
policy1(1:(end-1),1,:) = 1;

[value_function1,~,~,~] = SolvePDE3D(policy1);

policy2 = InitMBMPolicy();
policy2(1:(end-1),1,:) = 1;
policy2(2:end,3:end,:) = 2;

[value_function2,~,~,~] = SolvePDE3D(policy2);