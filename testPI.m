InitPara()

initial_policy = InitMBMPolicy();

[~,policy_PI,n_iter_PI] = PolicyIter3D(initial_policy);

[~,policy_VI,n_iter_VI] = ValueIter3D(policy_PI{n_iter_PI});

isequal(policy_PI{n_iter_PI},policy_VI)

sum(policy_PI{n_iter_PI}-policy_VI)