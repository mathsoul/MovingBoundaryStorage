InitPara()
global A_hold b_hold A_buy b_buy A_sell b_sell
[A_hold,b_hold] = GenerateHoldEquation();
[A_buy,b_buy] = GenerateBuyEquation();
[A_sell,b_sell] = GenerateSellEquation();

global NumQ NumX NumS

load('PI_policy_seq.mat')
load('MBM_policy_seq.mat')
load('MBM_sim_policy_seq.mat')

n_PI = length(policy_PI);
n_MBM = length(policy_MBM);
n_MBM_sim = length(policy_MBM_sim);


optimal_value = zeros(NumQ,NumX,NumS);
for i = 1:NumQ
    for j = 1:NumX
        for k = 1:NumS
            optimal_value(i,j,k) = max([value_function_PI{n_PI}(i,j,k),...
    value_function_MBM{n_MBM}(i,j,k),...
    value_function_MBM_sim{n_MBM_sim}(i,j,k)]);
        end
    end
end


diff_V_initial = log(norm(reshape(value_function_PI{1}-optimal_value,[NumQ*NumX*NumS,1]),Inf));

diff_V_PI = zeros(1,n_PI-1);

for i = 1:(n_PI-1)
    diff_V_PI(i) = log(norm(reshape(value_function_PI{i+1}-optimal_value,[NumQ*NumX*NumS,1]),Inf));
end

diff_V_MBM = zeros(1,n_MBM);

for i = 1:n_MBM
    diff_V_MBM(i) = log(norm(reshape(value_function_MBM{i}-optimal_value,[NumQ*NumX*NumS,1]),Inf));
end

diff_V_MBM_sim = zeros(1,n_MBM_sim);

for i = 1:n_MBM_sim
    diff_V_MBM_sim(i) = log(norm(reshape(value_function_MBM_sim{i}-optimal_value,[NumQ*NumX*NumS,1]),Inf));
end

% Plot value difference vs # iteration
fig_conv_rate = figure;
plot(linspace(1,(n_MBM+2)/2,n_MBM+1),[diff_V_initial,diff_V_MBM],'k',...
    1:n_MBM_sim+1,[diff_V_initial,diff_V_MBM_sim],'r',...
    1:n_PI,[diff_V_initial,diff_V_PI],'b')
title('log inf norm')
legend('MBM Alt','MBM Sim','PI')
print(fig_conv_rate,'MBM Alt vs MBM Sim vs PI','-dpng')
close all 

% How MBM moves

clear m

for i = 1:n_MBM
    figure
    BoundaryLimit(policy_MBM{i}(:,:,1))
    title(i)
    m(i) = getframe;
end

movie2avi(m,'MBM_3D.avi','compression','none','fps',3)

close all

% How MBM moves

clear m

for i = 1:n_MBM_sim
    figure
    BoundaryLimit(policy_MBM_sim{i}(:,:,1))
    title(i)
    m(i) = getframe;
end

movie2avi(m,'MBM_sim_3D.avi','compression','none','fps',3)

close all

% how PI moves
clear m

for i = 1:n_PI
    figure
    BoundaryLimit(policy_PI{i}(:,:,1))
    title(i)
    m(i) = getframe;
end

movie2avi(m,'PI_3D.avi','compression','none','fps',3)

close all
