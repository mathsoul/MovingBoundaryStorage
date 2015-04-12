% Function InitEquation() initialize the linear equations in the form that
% A x = b
% Inputs: None
% Outputs: 
%           A: matrix
%           b: vector


function [A,b] = InitEquation()
    global NumQ NumX NumS
    
    A = eye(NumQ,NumX,NumS);
    b = zeros(NumQ,NumX,NumS,1);
end