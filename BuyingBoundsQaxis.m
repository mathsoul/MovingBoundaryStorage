function buyingBound = BuyingBoundsQaxis(Policy)
NumX = length(Policy(1,:));
% NumQ = length(Policy(:,1));
buyingBound = zeros(NumX,1);


for i = 1:NumX
    if isempty(find(Policy(:,i) == 1,1,'last')) ~= 1
        buyingBound(i) = find(Policy(:,i) == 1,1,'last');
    end
end