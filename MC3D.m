InitPara()
global NumQ NumX NumS Smin Smax

[value_function_VI,optimal_policy_VI,n_iter_VI] = ValueIterwSeasonality();

disp_optimal_policy_VI = reshape4disp(optimal_policy_VI);

[value_function_PI,optimal_policy_PI,n_iter_PI] = PolicyIterwSeasonality(optimal_policy_VI);

disp_optimal_policy_PI = reshape4disp(optimal_policy_PI);


policy = reshape(optimal_policy_PI,NumQ,NumX,NumS);
% save('policySeasonality.mat','policy');
% load('policySeasonality.mat')

ds = (Smax - Smin)/(NumS - 1);
SVec = Smin:ds:Smax;

for i = NumS:-1:1
    s = SVec(i);
    
    figure
    BoundaryLimit(policy(:,:,i))
    title(round(s*2,2))
end
