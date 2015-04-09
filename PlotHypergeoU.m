a0 = beta/(2*kappa);
b0 = 1/2;

dq=(Qmax-Qmin)/(NumQ-1);       %Deltas
dx=(Xmax-Xmin)/(NumX-1);
XVec=Xmin:dx:Xmax;
QVec=Qmin:dq:Qmax;

% Nstart = find(XVec>alpha,1);
% 
% for i = 20:-1:1
%     U(21-i) = hypergeom(a0,b0,kappa/sigma^2*(XVec(i)-alpha)^2);
% end
% 
% for i = 2:20
%     R(i-1) = U(i)/U(i-1);
% end


% plot(1:20,U);
% 
% figure
% 
% plot(1:19,R);


for i = [1:NumX]
    U(i) = mchgu(a0,b0,kappa/sigma^2*(XVec(i)-alpha)^2);
    Uprime(i) = Hermite(-beta/(2*kappa),sqrt(kappa/sigma^2*(XVec(i)-alpha)^2));
end

