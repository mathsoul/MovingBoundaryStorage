% plotBoundary plots the boundary described in the buyingBoundary and sellingBoundary 
% function. Both functions are not under a wrong assumption. Therefore this
% function is useless.

function plotBoundary(alpha,beta,kappa,sigma,Qmax,Qmin,Xmax,Xmin)
    N = round((Xmax-Xmin)/0.01)+1;
    XVec = Xmin:0.01:Xmax;
    QBuyingVec = zeros(N,1);
    QSellingVec = zeros(N,1);
    for i = 1:N
        x = Xmin + (i-1) * 0.01;
        QBuyingVec(i) = buyingBoundary(x,alpha,beta,kappa,sigma,Qmax);
        QSellingVec(i) = sellingBoundary(x,alpha,beta,kappa,sigma,Qmax);
    end
    
    
    plot(QBuyingVec,XVec);
    
    hold on
    
    plot(Qmin:Qmax,alpha)
    
    plot(QSellingVec,XVec);
    

    
    hold off
     
end