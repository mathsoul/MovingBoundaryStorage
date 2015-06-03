InitPara()
global NumQ NumX NumS

load('PI_policy_seq.mat')
load('MBM_policy_seq0.1.mat')

n_PI = length(policy_PI);
n_MBM = length(policy_MBM);

optimal_policy = policy_PI{n_PI};
optimal_value =  SolvePDE3D(optimal_policy);

diff_V_initial = log(norm(reshape(value_function_PI{1}-optimal_value,[NumQ*NumX*NumS,1]),Inf));

diff_V_PI = zeros(1,n_PI-1);

for i = 1:(n_PI-1)
    diff_V_PI(i) = log(norm(reshape(value_function_PI{i+1}-optimal_value,[NumQ*NumX*NumS,1]),Inf));
end

diff_V_MBM = zeros(1,n_MBM);

for i = 1:n_MBM
    diff_V_MBM(i) = log(norm(reshape(value_function_MBM{i}-optimal_value,[NumQ*NumX*NumS,1]),Inf));
end

plot(1:(n_MBM+1),[diff_V_initial,diff_V_MBM],'k')
hold on
plot(1:n_PI,[diff_V_initial,diff_V_PI],'b')
hold off
title('log inf norm')
