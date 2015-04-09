function plotTransactionTimes(NPlot,BuyingTimes,SellingTimes,Xmax,Xmin,Qmax,Qmin,T,tau,N)
    NumX = length(BuyingTimes(:,1));
    NumQ = length(BuyingTimes(1,:));
    XVec = linspace(Xmin,Xmax,NumX);
        q = Qmin + (NPlot-1)*(Qmax-Qmin)/(NumQ-1);
%         figure
        plot(XVec,flipud(BuyingTimes(:,NPlot)),'--')
        hold on 
        plot(XVec,flipud(SellingTimes(:,NPlot)))
        hold off
        title({['Buying times and selling times'] ;['q =' num2str(q) ', dashed line is buying times, time limit is ' num2str(T)]; ['the length of each step is ' num2str(tau)];['the number of simulation is '  num2str(N)]})
        figure
        plot(XVec,flipud(SellingTimes(:,NPlot)+BuyingTimes(:,NPlot)))
        title({['Total transaction times'];['q =' num2str(q) ', time limit is ' num2str(T)]; ['the length of each step is ' num2str(tau)];['the number of simulation is '  num2str(N)]})
end