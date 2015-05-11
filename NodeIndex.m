% function NodeIndex() creates the node index in a matrix
% inputs: 
%           global variables NumQ (NumX, NumS) is the number of discretization on Q (X, S) axis
% outputs:
%           G is the node index matrix
%           Gi (Gj, Gk) is the q (x, s) index matrix

function [indexMat,indexVecQ,indexVecX,indexVecS] = NodeIndex()
    
    global NumQ NumX NumS 

    indexMat = reshape(1:NumQ*NumX*NumS,NumQ,NumX,NumS); % node index matrix

    indexVecQ = zeros(NumQ*NumX*NumS,1); 
    indexVecX = zeros(NumQ*NumX*NumS,1); 
    indexVecS = zeros(NumQ*NumX*NumS,1); 

    for l = 1:NumQ*NumX*NumS
        indexVecQ(l) = rem(l,NumQ);
        indexVecX(l) = ceil(rem(l,NumQ*NumX)/NumQ);
        indexVecS(l) = ceil(l/(NumQ * NumX));
    end

    indexVecQ(indexVecQ == 0) = NumQ;
    indexVecX(indexVecX == 0) = NumX;
end

