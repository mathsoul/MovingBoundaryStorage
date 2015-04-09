function plotTransactionCosts(NPlot,TotalBuyingCosts,TotalSellingCosts,Xmax,Xmin,Qmax,Qmin,T,tau,N)
    NumX = length(TotalBuyingCosts(:,1));
    NumQ = length(TotalBuyingCosts(1,:));
    XVec = linspace(Xmin,Xmax,NumX);
        q = Qmin + (NPlot-1)*(Qmax-Qmin)/(NumQ-1);
%         figure
        plot(XVec,flipud(TotalBuyingCosts(:,NPlot)),'--')
        hold on 
        plot(XVec,flipud(TotalSellingCosts(:,NPlot)))
        hold off
        title({['Buying costs and selling costs'];['q =' num2str(q) ', dashed line is total costs, time limit is ' num2str(T)]; ['the length of each step is ' num2str(tau)];['the number of simulation is '  num2str(N)]})
        figure
        plot(XVec,flipud(TotalBuyingCosts(:,NPlot))+TotalSellingCosts(:,NPlot),'--')
        title({['Total transaction costs'];['q =' num2str(q) 'time limit is ' num2str(T)]; ['the length of each step is ' num2str(tau)];['the number of simulation is '  num2str(N)]})
end