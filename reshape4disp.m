% function reshape4disp() reshape the vector for display
% inputs:
%       vec is the vector that is going to be reshaped
%       NumQ (NumX, NumS) is the number of discretization on Q (X, S) axis
% outputs:
%       array is the array after being reshaped

function array = Reshape4Disp(vec)
    global NumQ NumX NumS
    array = reshape(vec,NumQ,NumX,NumS);
    for k = 1:NumS
        array(:,:,k) = flipud(array(:,:,k)');
    end
end