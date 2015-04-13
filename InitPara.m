% function InitPara() initialize all the parameters as global variables 
% 
% inputs: None
% outputs: None


function InitPara()
    global kappa sigma alpha Qmax Qmin Xmax Xmin Smax Smin beta NumX NumQ...
        NumS BuyingType BuyingCostPara SellingType SellingCostPara ErrorTol ...
        DiffLevel MaxIteration
    
    kappa = 3.4;
    sigma = 0.59;
    alpha = 0.803;

    Qmax = 100;
    Qmin = 0;

    Xmax = 2.2;
    Xmin = -4;

    Smax = 1;
    Smin = -1;

    beta = 0.5;

    NumX = 21;
    NumQ = 21;
    NumS = 21;
    
    BuyingType = 'reciprocal';
    BuyingCostPara = [0.2,0,-0.2];

    SellingType = 'reciprocal';
    SellingCostPara = [0.2,0,-0.2];
    
    ErrorTol = 1;
    DiffLevel = 1e-3;
    MaxIteration = 1000;
end