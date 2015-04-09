%function plotValueFunction(NPlot,VPDE,VSim,Xmax,Xmin,Qmax,Qmin,T,tau,N)
%plots the value function from solving PDE and Simulation.


function plotValueFunction(NPlot,VPDE,VSim,StandardDeviation,Xmax,Xmin,Qmax,Qmin,T,tau,N)
    NumX = length(VPDE(:,1));
    NumQ = length(VPDE(1,:));
    XVec = linspace(Xmin,Xmax,NumX);
    UpperBound = flipud(VSim(:,NPlot) + 2*StandardDeviation(:,NPlot));
    LowerBound = flipud(VSim(:,NPlot) - 2*StandardDeviation(:,NPlot));
    q = Qmin + (NPlot-1)*(Qmax-Qmin)/(NumQ-1);



    plot(XVec,flipud(VPDE(:,NPlot)),'--')
    hold on 
    plot(XVec,flipud(VSim(:,NPlot)))
    for k = 1:NumX
    plot([XVec(k),XVec(k)],[UpperBound(k),LowerBound(k)],'-.')
    end
    hold off
    title({['q =' num2str(q) ', dashed line is from simulation, time limit is ' num2str(T)]; ['the length of each step is ' num2str(tau)];['the number of simulation is '  num2str(N)]})
end