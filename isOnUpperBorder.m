% Function isOnUpperBorder() test whether a node is on the upper border of
% a particular dimension
% Inputs:
%           index: node's index on a particular dimension
%           dimension_name
% Outputs:
%           indicator: whether this node is interior


function indicator = isOnUpperBorder(index,dimension_name)
    global NumX NumQ NumS
    
    if(dimension_name == 'Q')
        indicator = (index == NumQ);
    elseif(dimension_name == 'X')
        indicator = (index == NumX);
    else
        indicator = (index == NumS);
    end
end