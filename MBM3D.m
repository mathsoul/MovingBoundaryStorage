InitPara()

%% Move the boundary along X at the same time
initial_policy = InitMBMPolicy();

old_policy = initial_policy;

[value_function{1},profit_hold{1},profit_buy{1},profit_sell{1}] = SolvePDE3D(old_policy);  

new_policy{1} = MoveBoundaryAlongX3D(old_policy,profit_buy{1},profit_sell{1},'sell');

counter = 1;

while norm(new_policy{counter}-old_policy)>0
   old_policy = new_policy{counter};
   
   % Move buying boundary
   counter = counter +1;
   [value_function{counter},profit_hold{counter},profit_buy{counter},...
        profit_sell{counter}] = SolvePDE3D(old_policy);
   new_policy{counter} = ...
       MoveBoundaryAlongX3D(old_policy,profit_buy{counter},profit_sell{counter},'buy');
    
   % Move selling boundary
   counter = counter +1;
   [value_function{counter},profit_hold{counter},profit_buy{counter},...
        profit_sell{counter}] = SolvePDE3D(old_policy);
   new_policy{counter} = ...
       MoveBoundaryAlongX3D(old_policy,profit_buy{counter},profit_sell{counter},'sell');
end


MBM_policy = new_policy{counter};   
