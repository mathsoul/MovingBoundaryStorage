function y=Hermite(n,z)


%Slow
%y=2.^n.*hypergeomU(-n/2,1/2,z.^2);

%Fast
for i=1:length(z)
    if(z(i)<0)
y(i)=2^n*sqrt(pi)*(Hyper1F1(-n/2, 1/2, z.^2)/Gamma((1 - n)/2) - (2.*z.*Hyper1F1((1 - n)/2, 3/2, z.^2)/Gamma(-n/2)));
    else
       y(i)=2^n*mchgu(-n/2, 1/2, z(i)^2);
    end
end
