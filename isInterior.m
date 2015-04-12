% Function isInterior() test whether a node is interior in a particular
% dimension
% Inputs:
%           index: node's index on a particular dimension
%           dimension_name
% Outputs:
%           indicator: whether this node is interior


function  indicator = isInterior(index,dimension_name)
    global NumX NumQ NumS
    
    if(dimension_name == 'Q')
        indicator = (ismember(index, 2:(NumQ-1)));
    elseif(dimension_name == 'X')
        indicator = (ismember(index, 2:(NumX-1)));
    else
        indicator = (ismember(index, 2:(NumS-1)));
    end
end