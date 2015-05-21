load('VI.mat')
load('PI.mat')
load('MBM_0.1.mat')
% load('MBM.mat')
InitPara()
global NumS Smin Smax

ds = (Smax - Smin)/(NumS - 1);
SVec = Smin:ds:Smax;

for i = NumS:-1:1
    s = SVec(i);
    
    figure
    BoundaryLimit(policy_MBM(:,:,i))
    title({'MBM',round(s*2,2)})

    
    figure
    BoundaryLimit(policy_VI(:,:,i))
    title({'VI',round(s*2,2)})
    
    figure
    BoundaryLimit(policy_PI(:,:,i))
    title({'PI',round(s*2,2)})
end

mean(mean(mean(value_function_MBM-value_function_PI)))
mean(mean(mean(value_function_MBM-value_function_PI)))/mean(mean(mean(value_function_PI)))
