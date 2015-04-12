% Function isOnLowerBorder() test whether a node is on the upper border of
% a particular dimension
% Inputs:
%           index: node's index on a particular dimension
%           dimension_name
% Outputs:
%           indicator: whether this node is interior


function indicator = isOnLowerBorder(index,dimension_name)
    
    indicator = (index == 1);

end