InitPara()

[value_function_VI,optimal_policy_VI,n_iter_VI] = ValueIterwSeasonality();

Reshape4Disp(optimal_policy_VI);

[value_function_PI,optimal_policy_PI,n_iter_PI] = PolicyIterwSeasonality(optimal_policy_VI);

Reshape4Disp(optimal_policy_PI);