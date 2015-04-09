%This method is right when z is small but wrong when z is large. Actually,
%when z approaches infnity, it goes to infinity as well which contradicts
%the property hypergeomU function should have.

function [ v ] = hypergeomU(a,b,z);
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

v=pi/sin(pi*b)*(hypergeom(a,b,z)./(gamma(1+a-b)*gamma(b)) - z.^(1-b).*hypergeom(1+a-b,2-b,z)./(gamma(a)*gamma(2-b)));


%v=gamma(b-1)/gamma(a).*z.^(1-b).*hypergeom(a-b+1,2-b,z);
%v=v+gamma(1-b)/gamma(a-b+1).*hypergeom(a,b,z);



end

