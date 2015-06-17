% function InitPara() initialize all the parameters as global variables 
% 
% inputs: None
% outputs: None


function InitPara()
    global kappa sigma alpha Qmax Qmin Xmax Xmin Smax Smin Sratio beta NumX NumQ...
        NumS BuyingType BuyingCostPara SellingType SellingCostPara ErrorTol ...
        DiffLevel MaxIteration SellAbsDiff BuyAbsDiff CostPriceRatio
    
    kappa = 3.4;
    sigma = 0.59;
    alpha = 0.803;

    Qmax = 100;
    Qmin = 0;

    Xmax = 1.5;
    Xmin = 0;

    Smax = 5/4;
    Smin = 1/4;
    
    Sratio = 0.1;
    
    beta = 0.5;

    
    NumQ = 21;
    NumX = 21;
    NumS = 13;
    
    Smax = Smax - (Smax-Smin)/(NumS-1);
    NumS = NumS - 1;
    
    BuyingType = 'reciprocal';
    BuyingCostPara = [0.2,0,-0.2];

    SellingType = 'reciprocal';
    SellingCostPara = [0.2,0,-0.2];
    
    ErrorTol = 1;
    DiffLevel = 0;
    MaxIteration = 1000;
    
    SellAbsDiff = 0.0;
    BuyAbsDiff = 0.0;
    
    CostPriceRatio = 0.1;
end