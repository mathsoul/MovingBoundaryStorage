I = [0,2,2,2;0,0,0,0;0,1,1,0;1,1,1,0];
I = (flipud(I))';
[V2,Ih,Ib,Is,A,Aholding,Abuying,Aselling,C,Cholding,Cbuying,Cselling] = SolvePDE(alpha,beta,kappa,sigma,Xmax,Xmin,Qmax,Qmin,length(I(1,:)),length(I(:,1)),I);

for k =1: length(I(1,:))*length(I(:,1))
    inv(A(1:k,1:k))
end