% [Valuefunction, RegionIndicater, HoldingProfit, BuyingProfit,
% SellingProfit] =
% SolvePDE(kappa,sigma,alpha,beta,Xmax,Xmin,Qmax,Qmin,NumQ,NumX,lambda,mu)
% Computes the value function via solving PDE numerically by discretizing
% it. Here the boundary is fixed and analytical which is discribed in
% another two functions buyingBoundary and sellingBoundary.
%
% The inputs are:
%   alpha is the long term mean value
%   beta is the discount factor
%   kappa is the rate the value reverts to the mean
%   sigma is the volatility 
%   Xmax is the maximum value of log price
%   Xmin is the minimum value of log price
%   Qmax is the maximum value of volume
%   Qmin is the minimum value of volume
%   NumX is the number of pieces used to discretize log price
%   NumQ is the number of pieces used to discretize volume
%
% and returns:
%   Valuefunction is a NumX by NumQ matrix. Each component is the value of
%   begining with related log price and volume.
%   RegionIndicater is a NumX by NumQ matrix where each element tells us   
%   which regiton the related point belongs to. By using 1 representing 
%   buying, 2 representing selling and 0 representing holding.
%   HoldingProfit is a NumX by NumQ matrix where each element tells us   
%   what the instant profit is if holding at the related point.
%   BuyingProfit is a NumX by NumQ matrix where each element tells us   
%   what the instant profit is if buying at the related point.
%   SellingProfit is a NumX by NumQ matrix where each element tells us   
%   what the instant profit is if selling at the related point.
function RegionIndicater = ...
    RegionIndicater(alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,NumX,NumQ)
%% Create the grid

dq=(Qmax-Qmin)/(NumQ-1);       %Deltas
dx=(Xmax-Xmin)/(NumX-1);

%% Create Node index

Gi=reshape(repmat(1:NumQ,1,NumX),NumQ*NumX,1);                   % q index
Gj=reshape(repmat(1:NumX,NumQ,1),NumQ*NumX,1);                    % x index

RegionIndicater = zeros(NumQ,NumX);

    

for ij=1:NumQ*NumX    
    i=Gi(ij);j=Gj(ij);          %%Get i and j index
    q=(i-1)*dq+Qmin;
    x=(j-1)*dx+Xmin;
    
    if(q < buyingBoundary(x,alpha,beta,kappa,sigma,Qmax))
       RegionIndicater(i,j) = 1;
    elseif(q > sellingBoundary(x,alpha,beta,kappa,sigma,Qmax))
        RegionIndicater(i,j)=2;
    else
        RegionIndicater(i,j)=0;
    end

end

end


% % Checks whether the dimension of lambda and mu are both NumQ
% function[] = checkDimensions(lambda, mu, NumQ)
% 
% if length(lambda) ~= NumQ
%     error('lambda does not have the same dimension as NumQ')
% end
% 
% if length(mu) ~= NumQ
%     error('mu does not have the same dimension as NumQ')
% end
% end











