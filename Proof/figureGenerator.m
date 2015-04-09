% % This way generates terrible figure
% syms t
% x = 100 + 100*sin(t);
% y = 4 + 2*cos(t);
% 
% ezplot(x,y)
% xlim([0,100]);
% ylim([-4,4]);

n = 1000;
t = linspace(pi,3*pi/2,n);
alpha = 0.8;


%% Figure that shows the S^n, B^n H^n

figure

hold on 

xsellto = 100 + 100*sin(t);
ysellto = 2.3 + 1.2*cos(t);
 
xsellfrom = 100 + 100*sin(t);
ysellfrom = 2.3 + 0.6*cos(t);


plot(xsellto,ysellto,'-',xsellfrom,ysellfrom,'--','LineWidth',2)
xlim([0,100]);
ylim([-4,2.3]);

% t = linspace(0,2*pi,n);

% xbuyfrom = 0+ 60*sin(t);
% ybuyfrom = -2 + cos(t);

% plot(xbuyfrom,ybuyfrom)

plot([0,100],[alpha,alpha],'LineWidth',2)

gtext('S^n')
gtext('S^{n+1}')
% gtext('B^n = B^{n+1}')


% title('q = 0 is not isolated')
hold off

%% Figure that shows the delta V in stage 1

figure

hold on 

plot(xsellto,ysellto,'-',xsellfrom,ysellfrom,'--','LineWidth',2)
xlim([0,100]);
ylim([-4,2.3]);

% t = linspace(0,2*pi,n);

% xbuyfrom = 0+ 60*sin(t);
% ybuyfrom = -2 + cos(t);
% 
% plot(xbuyfrom,ybuyfrom)

plot([0,100],[alpha,alpha],'LineWidth',2)

gtext('\Delta V_q^{n+1}(x,q) = 0')
gtext('\Delta V_q^{n+1}(x,q) > 0')
% gtext('\Delta V_q^{n+1}(x,q) = 0')
gtext('L\Delta V^{n+1}(x,q) = 0')
% gtext('L\Delta V_q^{n+1}(x,q) = 0')

hold off

%% Figure that shows the delta V in stage 1

figure

hold on 

plot(xsellfrom,ysellfrom,'-')
xlim([0,100]);
ylim([-4,4]);

t = linspace(0,2*pi,n);

xbuyfrom = 0+ 60*sin(t);
ybuyfrom = -2 + cos(t);

xbuyto = 0 + 80*sin(t);
ybuyto = -2 + 1.5*cos(t);

plot(xbuyto,ybuyto,'-',xbuyfrom,ybuyfrom,'--')

plot([0,100],[alpha,alpha])

gtext('\Delta V_q^{n+1}(x,q) = 0')
gtext('\Delta V_q^{n+1}(x,q) < 0')
gtext('\Delta V_q^{n+1}(x,q) = 0')
gtext('L\Delta V^{n+1}(x,q) = 0')


hold off

%% Situation 1 where q = 0 is isolated state


xsellfrom = 100 + 100*sin(t);
ysellfrom = 4 + 2*cos(t);

xsellto = 100 + 100*sin(t);
ysellto = 4 + 3*cos(t);

hold on 

plot(xsellto,ysellto,'-',xsellfrom,ysellfrom,'--')
xlim([0,100]);
ylim([-4,4]);

t = linspace(0,2*pi,n);

xbuyfrom = 50 + 30*sin(t);
ybuyfrom = -2 + cos(t);

plot(xbuyfrom,ybuyfrom)

plot([0,100],[alpha,alpha])

% title('q = 0 is isolated')

hold off

%% Situation 2 where q = 0 is not isolated

figure

hold on 

plot(xsellto,ysellto,'-',xsellfrom,ysellfrom,'--')
xlim([0,100]);
ylim([-4,4]);

t = linspace(0,2*pi,n);

xbuyfrom = 0+ 60*sin(t);
ybuyfrom = -2 + cos(t);

plot(xbuyfrom,ybuyfrom)

plot([0,100],[alpha,alpha])

% title('q = 0 is not isolated')
hold off