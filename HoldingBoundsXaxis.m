function [upperBound,lowerBound] = HoldingBoundsXaxis(Policy)
NumX = length(Policy(1,:));
NumQ = length(Policy(:,1));
upperBound = zeros(NumQ,1) + NumX + 1;
lowerBound = zeros(NumQ,1);

for i = 1:NumQ
    if isempty(find(Policy(i,:) == 2,1,'first')) ~= 1
    	upperBound(i) = find(Policy(i,:) == 2,1,'first');
    end
    if isempty(find(Policy(i,:) == 1,1,'last')) ~=1
        lowerBound(i) = find(Policy(i,:) == 1,1,'last');
    end
end