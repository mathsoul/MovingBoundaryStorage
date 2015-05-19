load('VI.mat')
load('PI.mat')
% load('MBM_0.1.mat')
% load('MBM.mat')
global NumS Smin Smax

ds = (Smax - Smin)/(NumS - 1);
SVec = Smin:ds:Smax;

for i = NumS:-1:1
    s = SVec(i);
    
%     figure
%     BoundaryLimit(MBM_policy(:,:,i))
%     title({'MBM',round(s*2,2)})
    
    figure
    BoundaryLimit(optimal_policy_VI(:,:,i))
    title({'VI',round(s*2,2)})
    
    figure
    BoundaryLimit(optimal_policy_PI(:,:,i))
    title({'PI',round(s*2,2)})
end

% mean(mean(mean(value_function_MBM-value_function_PI)))
% mean(mean(mean(value_function_MBM-value_function_PI)))/mean(mean(mean(value_function_PI)))