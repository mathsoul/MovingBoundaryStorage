
%% Parameters

% sig=.334;  %vol
% kappa=.301;  % rev speed
% gamma=3.093; %rev level for (ln S)
% lambda= 0; %mkt price of risk
% al=(gamma - sig^2/(2*kappa)) - lambda;  %risk neutral rev level for X
% beta=.3;
% 
% 
% Smid=exp(gamma);
% Smax=50;
% Smin=.1;
% Xmax=log(Smax);
% Xmin=-3; %log(Smin);


kappa = 3.4;
sig =0.59;
al = 0.802;
beta = 0.05;

% syms x

Xmax= 3  ;
% Xmin=double(solve(1/2*sigma^2*exp(x)+kappa*(alpha-x-1)*exp(x) - beta*(exp(x))==0,x))-1;
Xmin = -1.5;


Qmax=100;
Qmin=0;

Ni = 21;
Nj = 21;


c1=beta/(2*kappa);
c2=0.5;
c3=kappa/sig^2;



%% Grid Parameters

% Ni=26;          %q disc
% Nj=Ni;          %x disc





%% costs
% lambda=(1.05).^(([0:Ni-1]))/1;
% mu=fliplr(lambda)*10;
% mu=ones(1,Ni)*16;



%%output filename
filename=['Storage' num2str(Ni) '-' num2str(Nj) ];





%% Grid pre-calcs

dq=(Qmax-Qmin)/(Ni-1);       %Deltas
dx=(Xmax-Xmin)/(Nj-1);
XVec=Xmin:dx:Xmax;
QVec=Qmin:dq:Qmax;


lambda = buyingCost(QVec,Qmax);
mu = sellingCost(QVec,Qmax);

G=repmat(0:Ni:Ni*Nj-1,Ni,1)+repmat(1:Ni,Nj,1)';        %Node index mat
Gi=reshape(repmat(1:Ni,1,Nj),Ni*Nj,1);                   % q index
Gj=reshape(repmat(1:Nj,Ni,1),Ni*Nj,1);                    % x index




%%


%% Setting init values


%b(1,:)=[ones(1,Ni-1)*1 0];  %% Guess Boundaries
%s(1,:)=[Nj+1 ones(1,Ni-1)*Nj];

%b(1,:)=[ones(1,Ni)*1 ];  %% Guess Boundaries
%s(1,:)=[ones(1,Ni)*Nj];

%b0=4;b1=1;
%s0=Nj;s1=Nj-3;

%b(1,:)=round([b0:-(b0-b1)/(Ni-1):b1]);
%s(1,:)=round([s0:-(s0-s1)/(Ni-1):s1]);


%%



%%
%%%%%%%%%%%%%%%%%% Computing begins %%%%%%%%%%%%%%%

tic
disp('----------------------------------------------------');


Bit=1;                              %% Boudnary iteration index
BitMax=10;                          %% Max iterations
Run=1;                              %% Convergence flag




%A=sparse(zeros(Ni*Nj));         %% Initialize
%Ab=sparse(zeros(Ni*Nj));         %% Initialize
%As=sparse(zeros(Ni*Nj));         %% Initialize

A=zeros(Ni*Nj);         %% Initialize
Ab=zeros(Ni*Nj);         %% Initialize
As=zeros(Ni*Nj);         %% Initialize

C=zeros(Ni*Nj,1);
Cb=zeros(Ni*Nj,1);
Cs=zeros(Ni*Nj,1);


% Assemble


for ij=1:Ni*Nj
    eq=ij;
    
    i=Gi(ij);j=Gj(ij);          %%Get i and j index
    q=(i-1)*dq+Qmin;
    x=(j-1)*dx+Xmin;
    
    
    Cxx=0.5*sig^2/dx^2;
    Cx=kappa*(al-x)/dx;
    
    
    %hold
    Cden=(2*Cxx+abs(Cx)+beta);
    if(j>1 && j<Nj)
        
        A(eq,G(i,j+1))=Cxx/Cden;
        A(eq,G(i,j-1))=Cxx/Cden;
        if(Cx>0)
            A(eq,G(i,j+1))=A(eq,G(i,j+1))+Cx/Cden;
        else
            A(eq,G(i,j-1))=A(eq,G(i,j-1))-Cx/Cden;
        end
        C(eq)=0;
        
    elseif(j==1)
        
        
        %A(eq,G(i,j+1))=1;
        
        %bp=max(0,Cx);bm=min(0,Cx);
        %A(eq,G(i,j+2))=Cxx/(beta+bp);
        %A(eq,G(i,j+1))=(-2*Cxx+Cx)/(beta+bp);
        %A(eq,G(i,j))=-bm; %most likely 0
        %if (-2*Cxx+Cx)<0
        %    disp('Warning: lower truncation probably needs to be brought down');
        %end
        
        
%         A(eq,G(i,j+1))=hypergeomU(c1,c2,c3*(XVec(j)-al)^2)/hypergeomU(c1,c2,c3*(XVec(j+1)-al)^2);
        A(eq,G(i,j+1))=Hermite(-2*c1,sqrt(c3)*abs(XVec(j)-al))/Hermite(-2*c1,sqrt(c3)*abs(XVec(j+1)-al));
        C(eq)=0;
        
    elseif(j==Nj)
        %A(eq,G(i,j-1))=1;
        
        
        
        %bp=max(0,Cx);bm=min(0,Cx);
        %A(eq,G(i,j-2))=Cxx/(beta-bm);
        %A(eq,G(i,j-1))=(-2*Cxx-Cx)/(beta-bm);
        %A(eq,G(i,j))=bp; %most likely 0
        %if (-2*Cxx-Cx)<0
        %    disp('Warning: upper truncation probably needs to be pullsed up');
        %end
        
%         A(eq,G(i,j-1))=hypergeomU(c1,c2,c3*(XVec(j)-al)^2)/hypergeomU(c1,c2,c3*(XVec(j-1)-al)^2);
        A(eq,G(i,j-1))=Hermite(-2*c1,sqrt(c3)*abs(XVec(j)-al))/Hermite(-2*c1,sqrt(c3)*abs(XVec(j-1)-al));
        C(eq)=0;
    end
    
    %buy
    if(i<Ni)
        Ab(eq,G(i+1,j))=1;
        Cb(eq)=-(exp(x)+lambda(i))*dq;
    end
    
    %sell
    if(i>1)
        As(eq,G(i-1,j))=1;
        Cs(eq)=(exp(x)-mu(i))*dq;
    end
    
end



if(0) %value iteration
    W{1}=zeros(Ni*Nj,1);
    run=1;k=1;
    while(run)
        [a b]=max([A*W{k}+C,Ab*W{k}+Cb,As*W{k}+Cs],[],2);
        W{k+1}=a; p{k+1}=b;
        disp(norm(W{k+1}-W{k}));
        if(norm(W{k+1}-W{k})<0.001)
            run=0;
        end
        k=k+1;
    end
else %policy iteration
    W{1}=zeros(Ni*Nj,1);p{1}=ones(Ni*Nj,1);
    run=1;k=1;
    while(run && k<10000)
        
        if(k>5)
            1
        end
        
        Cn=C.*(p{k}==1)+Cb.*(p{k}==2)+Cs.*(p{k}==3);
        An=A.*repmat((p{k}==1),1,Ni*Nj)+Ab.*repmat((p{k}==2),1,Ni*Nj)+As.*repmat((p{k}==3),1,Ni*Nj);
        W{k}=(eye(Ni*Nj)-An)\Cn;
        
        
        [a b]=max([A*W{k}+C,Ab*W{k}+Cb,As*W{k}+Cs],[],2);
        p{k+1}=b;
        
        
        disp(norm(p{k+1}-p{k}));
        if(norm(p{k+1}-p{k})<0.0001)
            run=0;
        end
        
        k=k+1;
    end
    
    k=k-1;
    
end



bl=ones(Nj,1);
bu=Ni*ones(Nj,1);


figure;
mesh((reshape(Gi,Ni,Nj)-1)*dq+Qmin,(reshape(Gj,Ni,Nj)-1)*dx+Xmin,reshape(W{k},Ni,Nj));
ax=axis;
ax(4)=dx*(Nj-1)+Xmin;
axis(ax);
xlabel('Q');
ylabel('X');

figure; hold on;
for j=1:Nj
    for i=1:Ni
        if(p{k}(G(i,j))==1)  %hold
            plot((i-1)*dq+Qmin,(j-1)*dx+Xmin,'r.');
        elseif(p{k}(G(i,j))==2)  %buy
            plot((i-1)*dq+Qmin,(j-1)*dx+Xmin,'b.');
            if(bl(j)<i)
                bl(j)=i;
            end
        elseif(p{k}(G(i,j))==3)  %sell
            plot((i-1)*dq+Qmin,(j-1)*dx+Xmin,'g.');
            if(bu(j)>i)
                bu(j)=i;
            end
        end
    end
end
line([Qmin Qmax],[al al]);

if(0)
    figure; hold;
    line((bl-1)*dq+Qmin,[0:Ni-1]'*dx+Xmin);
    line((bu-1)*dq+Qmin,[0:Ni-1]'*dx+Xmin);
    axis([Qmin Qmax Xmin Xmax]);
    line([Qmin Qmax Qmax Qmin],[Xmin Xmin Xmax Xmax]);
    xlabel('Q');
    ylabel('X');
end

