
function plotGrids(Xmax,Xmin,Qmax,Qmin,NumX,NumQ)
    X = linspace(Xmin,Xmax,NumX);
    Q = linspace(Qmin,Qmax,NumQ);
    
    QGrid1 = ones(NumX,1)*Q;
    XGrid1 = X'*ones(1,NumQ);
    
    QGrid2 = Q'*ones(1,NumX);
    XGrid2 = ones(NumQ,1)*X;
    
    hold on
    
    plot(QGrid1,XGrid1,'b')
    
    plot(QGrid2,XGrid2,'b')
    
    hold off
end