%logPriceGenerator(x,alpha,kappa,sigma,T,tau) generates the log price with
%initial price x.

function X = logPriceGenerator(x,alpha,kappa,sigma,tau,N)
%     x = 0;
%     tau = 0.01;
%     N = 10000;
    X = zeros(N,1); %log price vector w.r.t time
    X(1) = x;
    for k = 1:(N-1)
        Z = normrnd(0,1,1,1);
        X(k+1) = exp(-kappa*tau) * X(k) + alpha * ( 1 - exp(-kappa*tau)) + sqrt((1-exp(-2*kappa*tau))/(2*kappa))*sigma * Z;
%         if(X(k+1)> 1.8)
%             Num1 = k+1; break
%         elseif(X(k+1)>1.6 )
%             Num2 = k+1;
%         end
    end
%     disp(max(X))
end