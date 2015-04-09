[x,y] = meshgrid([-1.05:.2:3.75]);
z = x.*exp(-x.^2-y.^2);
axis tight;
set(gca,'nextplot','replacechildren');
for j = 1:40
      surf(x*sin(pi*j/100),y*sin(pi*j/100),z*sin(-pi*j/100));
      m(j) = getframe;
end
movie(m)