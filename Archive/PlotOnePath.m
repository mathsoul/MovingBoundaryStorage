%Function plotOnePath(x,q,alpha,beta,kappa,sigma,Qmax,Qmin,T,tau,N) gives
%the path with N points and all the input parameters.


function PlotOnePath(x,q,alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,T,tau,N)
%     figure
    
    plotBoundary(alpha,beta,kappa,sigma,Qmax,Qmin,Xmax,Xmin)
    
    hold on
    
    X = logPriceGenerator(x,alpha,kappa,sigma,tau,N); %log price vector w.r.t time
    Q = zeros(N,1); %volume vector w.r.t. time
    Q(1) = q;
        if( Q(1) <  buyingBoundary(X(1),alpha,beta,kappa,sigma,Qmax)) % buy at initial time
            Q(1) =buyingBoundary(X(1),alpha,beta,kappa,sigma,Qmax);
            scatter(q,X(1),5);%price before trading
            plot([q,Q(1)],[X(1),X(1)],'g');
        elseif(Q(1)> sellingBoundary(X(1),alpha,beta,kappa,sigma,Qmax)) %Sell at initial time
            Q(1) = sellingBoundary(X(1),alpha,beta,kappa,sigma,Qmax);
            scatter(q,X(1),5);%price before trading
            plot([q,Q(1)],[X(1),X(1)],'r');
        end
    for k = 1:(N-1)
%         Buy, sell or hold
        if( Q(k) <  buyingBoundary(X(k+1),alpha,beta,kappa,sigma,Qmax)) % buy
            Q(k+1) =buyingBoundary(X(k+1),alpha,beta,kappa,sigma,Qmax);
            scatter(Q(k),X(k+1),5);%price before trading
            plot([Q(k),Q(k+1)],[X(k+1),X(k+1)],'g');
        elseif(Q(k)> sellingBoundary(X(k+1),alpha,beta,kappa,sigma,Qmax)) %Sell
            Q(k+1) = sellingBoundary(X(k+1),alpha,beta,kappa,sigma,Qmax);
            scatter(Q(k),X(k+1),5);%price before trading
            plot([Q(k),Q(k+1)],[X(k+1),X(k+1)],'r');
        else
            Q(k+1) = Q(k);
        end
    end
    
    scatter(Q,X,5) %price after trading
    
    title(['Starting point is (' num2str(x) ',' num2str(q) ')'])
    
    hold off

end