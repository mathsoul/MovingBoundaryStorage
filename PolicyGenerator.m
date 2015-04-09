%Function PolicyGenerator(BuyingBoundsQ,SellingBoundsQ,NumX,NumQ) 
%generates the policy from Buying bounds along Q and Selling bounds along Q

function Policy = PolicyGenerator(BuyingBoundsQ,SellingBoundsQ,NumX,NumQ)
    
    Policy = zeros(NumQ,NumX);
    for j = 1 : NumX
        Policy(1:BuyingBoundsQ(j),j) = 1;
        Policy(SellingBoundsQ(j):NumQ,j) = 2;
    end
end