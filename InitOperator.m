% Function InitOperator() initialize the linear equations in the form that
% V(i,j,k) = operator V + constant
% Inputs: None
% Outputs: 
%           operator
%           constant


function [operator,constant] = InitOperator()
    global NumQ NumX NumS
    
    operator = zeros(NumQ*NumX*NumS);
    constant = zeros(NumQ*NumX*NumS,1);
end