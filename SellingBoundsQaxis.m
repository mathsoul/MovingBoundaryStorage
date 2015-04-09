function SellBound = SellingBoundsQaxis(Policy)
NumX = length(Policy(1,:));
NumQ = length(Policy(:,1));
SellBound = zeros(NumX,1) + NumQ + 1;

for i = 1:NumX
    if isempty(find((Policy(:,i) == 2),1,'first')) ~= 1
        SellBound(i) = find((Policy(:,i) == 2),1,'first');
    end
end